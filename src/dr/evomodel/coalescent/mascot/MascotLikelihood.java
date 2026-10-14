/*
 * MascotLikelihood.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.coalescent.AbstractStructuredCoalescentLikelihood;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentTipData;
import dr.evomodel.coalescent.StructuredTipStates;
import dr.evomodel.coalescent.basta.AbstractPopulationSizeModel;
import dr.evomodel.coalescent.basta.StructuredCoalescentLikelihoodGradient;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.ProcessAlongTree;
import dr.evomodel.treedatalikelihood.markovjumps.CompleteHistoryTraits;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

public class MascotLikelihood extends AbstractStructuredCoalescentLikelihood
        implements Citable, ProcessAlongTree, MarkovJumpsTraitProvider {

    public static final String MASCOT_LIKELIHOOD = "mascotLikelihood";

    /** Tree-trait tag for adjoint node-state sensitivities, not probabilities. */
    public static final String DEFAULT_ANCESTRAL_STATE_TAG_NAME = "mascot.stateSensitivity";

    public static final String DEFAULT_MODE_STATE_TAG_NAME = "mascot.state";

    /** Suffixes of the stochastic-mapping traits; the same names BASTA's Markov-jump add-ons use. */
    public static final String HISTORY_SUFFIX = ".history";
    public static final String COUNT_SUFFIX = ".count";
    public static final String RAW_SUFFIX = ".raw";

    private final TreeModel treeModel;
    private final PatternList tipPatterns;
    private final MascotDynamics dynamics;
    private final double maxStep;
    private final boolean checkProbabilities;
    private final String ancestralStateTagName;

    private MascotPreparedInput preparedInput;
    private MascotLikelihoodDelegate likelihoodDelegate;
    private boolean preparedInputKnown;
    private boolean likelihoodDelegateKnown;
    private boolean gradientKnown;
    private StructuredCoalescentTipData tipDataCache;

    private double[] thetaBuffer;
    private double[] branchRateBuffer;

    private boolean ancestralStatesKnown;

    // Stochastic mapping: the preparation (taped forward pass + Pi_h filter) lives in the delegate and is valid
    // until the model changes or is restored; the draw is also discarded whenever the state is stored, which the
    // MCMC does once per iteration right after logging, so every log step gets a fresh draw shared by all loggers.
    private boolean stochasticMapPrepared;
    private MascotStochasticMap stochasticMap;
    private boolean warnedStochasticMapLogLikelihood;
    private String stochasticMappingTag;
    private StochasticMapSampler stochasticMapSampler = StochasticMapSampler.GRID;
    private boolean warnedStiffFallback;

    private static final class DerivedBuffers {
        double[] combinedGradient;
        double[] clockGradient;
        double[] ancestralStateScores;
    }

    private final DerivedBuffers[] derivedBuffers = {new DerivedBuffers(), new DerivedBuffers()};
    private int currentDerivedIndex;
    private int storedDerivedIndex;

    private MascotPreparedInput storedPreparedInput;
    private MascotLikelihoodDelegate storedLikelihoodDelegate;
    private boolean storedPreparedInputKnown;
    private boolean storedLikelihoodDelegateKnown;
    private boolean storedGradientKnown;
    private boolean storedAncestralStatesKnown;

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            Parameter migrationRates,
                            SubstitutionModel[] migrationModels,
                            BranchRateModel branchRateModel,
                            AbstractPopulationSizeModel populationSizeModel,
                            int stateCount,
                            Parameter epochTimes,
                            double maxStep,
                            boolean checkProbabilities,
                            String ancestralStateTagName,
                            String modeStateTagName) {
        super(name == null ? MASCOT_LIKELIHOOD : name, treeModel, tipPatterns, stateCount, branchRateModel,
                migrationModels, populationSizeModel, tipPatterns.getDataType(), modeStateTagName);
        this.treeModel = treeModel;
        this.tipPatterns = tipPatterns;
        this.dynamics = migrationRates != null ?
                new MascotDynamics(stateCount, migrationRates, populationSizeModel, epochTimes) :
                new MascotDynamics(stateCount, migrationModels, populationSizeModel, epochTimes);
        this.maxStep = maxStep;
        this.checkProbabilities = checkProbabilities;
        this.ancestralStateTagName = ancestralStateTagName;
        if (migrationRates != null) {
            addVariable(migrationRates);
        } else {
            for (SubstitutionModel migrationModel : migrationModels) {
                addModel(migrationModel);
            }
        }
        addModel(populationSizeModel);
        if (epochTimes != null) {
            addVariable(epochTimes);
        }

        this.preparedInputKnown = false;
        this.likelihoodDelegateKnown = false;
        this.gradientKnown = false;

        treeTraits.addTrait(new TreeTrait.DA() {
            public String getTraitName() {
                return ancestralStateTagName;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public double[] getTrait(Tree tree, NodeRef node) {
                if (tree.isExternal(node)) {
                    return StructuredTipStates.getPartials(tree, node, tipPatterns, getStateCount(), true,
                            ancestralStateTagName + " tip trait");
                }
                double[] row = new double[getStateCount()];
                getAncestralStateScores(node.getNumber(), row);
                return row;
            }
        });

        addModeStateTrait(tag, this::getModeState);
    }

    public int getModeState(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            double[] tipPartials = StructuredTipStates.getPartials(tree, node, tipPatterns, getStateCount(), true,
                    tag + " tip trait");
            return argmax(tipPartials);
        }
        double[] row = new double[getStateCount()];
        getAncestralStateScores(node.getNumber(), row);
        return argmax(row);
    }

    public double[] getGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        return currentDerived().combinedGradient.clone();
    }

    public double[] getMigrationGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        double[] result = new double[dynamics.getMigrationRates().getDimension()];
        dynamics.writeMigrationGradient(currentDerived().combinedGradient, result);
        return result;
    }

    public double[] getPopSizeGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        double[] result = new double[dynamics.getPopSizes().getDimension()];
        dynamics.writePopSizeGradient(currentDerived().combinedGradient, result);
        return result;
    }

    public double[] getClockGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        double[] clockGradient = currentDerived().clockGradient;
        return clockGradient == null ? null : clockGradient.clone();
    }

    /** Node-major internal-node sensitivities, not posterior probabilities. */
    public double[] getAncestralStateScores() {
        if (!ancestralStatesKnown) {
            evaluateAncestralStates();
        }
        return currentDerived().ancestralStateScores.clone();
    }

    public void getAncestralStateScores(int nodeNumber, double[] destination) {
        if (destination.length != getStateCount()) {
            throw new IllegalArgumentException("destination length " + destination.length +
                    " does not match stateCount " + getStateCount());
        }
        if (!ancestralStatesKnown) {
            evaluateAncestralStates();
        }
        double[] ancestralStateScores = currentDerived().ancestralStateScores;
        if (nodeNumber < 0 || (nodeNumber + 1) * getStateCount() > ancestralStateScores.length) {
            throw new IllegalArgumentException("node number out of range: " + nodeNumber);
        }
        System.arraycopy(ancestralStateScores, nodeNumber * getStateCount(), destination, 0, getStateCount());
    }

    /**
     * One joint draw of every node's location and every branch's complete location history from Pi_h, MASCOT's
     * frozen-environment surrogate on its RK4 grid, given the current tree, tip data and parameter values (see
     * {@link MascotStochasticMapper}). Every caller gets the same draw until the model changes, the state is stored
     * or restored, or {@link #forceRedraw()} is called.
     *
     * @throws IllegalStateException if the log likelihood is not finite at the current state
     */
    public MascotStochasticMap getStochasticMap() {
        if (stochasticMap == null) {
            if (!stochasticMapPrepared) {
                prepareStochasticMap();
            }
            stochasticMap = likelihoodDelegate.drawStochasticMap();
            if (stochasticMap.getSampler() != stochasticMapSampler && !warnedStiffFallback) {
                Logger.getLogger("dr.evomodel").warning(String.format("MASCOT stochastic mapping: the RK4 steps are " +
                        "too stiff for the h-transform sampler (max h mu = %.2f > 1), so this state, and any other " +
                        "as stiff, is drawn with the exact grid sampler. A smaller maxStep lets the h-transform run.",
                        stochasticMap.getMaxRateTimesLength()));
                warnedStiffFallback = true;
            }
        }
        return stochasticMap;
    }

    /** A copy of the drawn history of the branch above {@code node}; null for the root. */
    public StateHistory getBranchHistory(NodeRef node) {
        return getStochasticMap().getBranchHistory(node.getNumber());
    }

    /** Discards the current draw, keeping its preparation, so the next request draws afresh. */
    public void forceRedraw() {
        stochasticMap = null;
    }

    /** Chooses how draws are made; {@link StochasticMapSampler#GRID} by default. */
    public void setStochasticMapSampler(StochasticMapSampler sampler) {
        if (sampler == null) {
            throw new IllegalArgumentException("sampler must not be null");
        }
        if (sampler != stochasticMapSampler) {
            stochasticMapSampler = sampler;
            stochasticMapPrepared = false;
            stochasticMap = null;
        }
    }

    public StochasticMapSampler getStochasticMapSampler() {
        return stochasticMapSampler;
    }

    /**
     * Registers the stochastic-mapping tree traits under {@code tag}, with BASTA's names and types: {@code tag}
     * (drawn node states), {@code tag.history} (each branch's changes as {site, height, from, to}), {@code tag.count}
     * (jumps per branch), {@code tag.raw} (copies of the {@link StateHistory} objects; not loggable) and
     * {@code tag.<register name>} per count (K x K) or reward (K) register; plus the untagged {@code history} and
     * {@code allTransitions} traits that completeHistoryLogger reads. Labels are backward in time: a change from a
     * to b (register entry a K + b) is a forward-time migration from b to a.
     */
    public void enableStochasticMapping(String tag, List<Parameter> countRegisters, List<Parameter> rewardRegisters) {
        if (stochasticMappingTag != null) {
            throw new IllegalStateException("stochastic mapping is already enabled under tag " + stochasticMappingTag);
        }
        if (tag == null || tag.isEmpty()) {
            throw new IllegalArgumentException("stochastic mapping needs a non-empty tag");
        }
        int K = getStateCount();
        List<String> names = new ArrayList<>();
        names.add(tag);
        names.add(tag + HISTORY_SUFFIX);
        names.add(tag + COUNT_SUFFIX);
        names.add(tag + RAW_SUFFIX);
        for (Parameter register : countRegisters) {
            if (register.getDimension() != K * K) {
                throw new IllegalArgumentException("count register " + register.getParameterName() + " has dimension " +
                        register.getDimension() + "; expected " + (K * K));
            }
            names.add(tag + "." + register.getParameterName());
        }
        for (Parameter register : rewardRegisters) {
            if (register.getDimension() != K) {
                throw new IllegalArgumentException("reward register " + register.getParameterName() + " has dimension " +
                        register.getDimension() + "; expected " + K);
            }
            names.add(tag + "." + register.getParameterName());
        }
        TreeTrait[] loggerTraits = CompleteHistoryTraits.createLoggerTraits(
                historyTrait(tag + HISTORY_SUFFIX), countTrait(tag + COUNT_SUFFIX));
        for (TreeTrait trait : loggerTraits) {
            names.add(trait.getTraitName());
        }
        Set<String> unique = new HashSet<>();
        for (String name : names) {
            if (!unique.add(name) || treeTraits.getTreeTrait(name) != null) {
                throw new IllegalArgumentException("stochastic-mapping trait name '" + name + "' is already in use");
            }
        }

        treeTraits.addTrait(new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public int[] getTrait(Tree tree, NodeRef node) {
                requireOwnTree(tree);
                return new int[]{getStochasticMap().getNodeState(node.getNumber())};
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formattedState(getTrait(tree, node));
            }
        });
        treeTraits.addTrait(historyTrait(tag + HISTORY_SUFFIX));
        treeTraits.addTrait(countTrait(tag + COUNT_SUFFIX));
        treeTraits.addTrait(new TreeTrait<StateHistory[]>() {
            public String getTraitName() {
                return tag + RAW_SUFFIX;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public Class getTraitClass() {
                return StateHistory[].class;
            }

            public StateHistory[] getTrait(Tree tree, NodeRef node) {
                StateHistory history = branchHistory(tree, node);
                return history == null ? null : new StateHistory[]{history};
            }

            public String getTraitString(Tree tree, NodeRef node) {
                throw new UnsupportedOperationException(tag + RAW_SUFFIX + " is not loggable");
            }

            public boolean getLoggable() {
                return false;
            }
        });
        for (Parameter register : countRegisters) {
            treeTraits.addTrait(registerTrait(tag + "." + register.getParameterName(), register.getParameterValues(),
                    true));
        }
        for (Parameter register : rewardRegisters) {
            treeTraits.addTrait(registerTrait(tag + "." + register.getParameterName(), register.getParameterValues(),
                    false));
        }
        treeTraits.addTraits(loggerTraits);
        stochasticMappingTag = tag;
        Logger.getLogger("dr.evomodel").info("MASCOT stochastic mapping: one joint draw of node states and branch " +
                "histories per log step, as traits '" + tag + "', '" + tag + HISTORY_SUFFIX + "', '" + tag +
                COUNT_SUFFIX + "' (" + stochasticMapSampler.getXmlName() + " sampler)");
    }

    private TreeTrait<String[]> historyTrait(final String name) {
        return new TreeTrait.SA() {
            public String getTraitName() {
                return name;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public String[] getTrait(Tree tree, NodeRef node) {
                StateHistory history = branchHistory(tree, node);
                return history == null ? null : new String[]{history.toStringChanges(1, dataType)};
            }
        };
    }

    private TreeTrait<double[]> countTrait(final String name) {
        return new TreeTrait.DA() {
            public String getTraitName() {
                return name;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public double[] getTrait(Tree tree, NodeRef node) {
                StateHistory history = branchHistory(tree, node);
                return history == null ? null : new double[]{history.getNumberOfJumps()};
            }
        };
    }

    /** Registered jumps (K x K register) or reward (K register) on each branch. */
    private TreeTrait<double[]> registerTrait(final String name, final double[] register, final boolean counts) {
        return new TreeTrait.DA() {
            public String getTraitName() {
                return name;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public double[] getTrait(Tree tree, NodeRef node) {
                StateHistory history = branchHistory(tree, node);
                if (history == null) {
                    return null;
                }
                return new double[]{counts ? history.getTotalRegisteredCounts(register) :
                        history.getTotalReward(register)};
            }
        };
    }

    private StateHistory branchHistory(Tree tree, NodeRef node) {
        requireOwnTree(tree);
        return tree.isRoot(node) ? null : getStochasticMap().getBranchHistory(node.getNumber());
    }

    private void requireOwnTree(Tree tree) {
        if (tree != treeModel) {
            throw new IllegalArgumentException("MASCOT stochastic-mapping traits are only defined on the likelihood's " +
                    "own tree");
        }
    }

    @Override
    public void calculatePostOrderStatistics() {
        getLogLikelihood();
    }

    @Override
    public String formattedState(int[] state) {
        return formattedState(state[0]);
    }

    private void prepareStochasticMap() {
        double logLikelihood = getLogLikelihood();
        if (!Double.isFinite(logLikelihood)) {
            throw new IllegalStateException("MASCOT stochastic map cannot be drawn: the log likelihood is " +
                    logLikelihood + " at the current state");
        }
        ensurePreparedInput();
        ensureLikelihoodDelegate();
        double[] nodeHeights = new double[treeModel.getNodeCount()];
        for (int i = 0; i < nodeHeights.length; i++) {
            NodeRef node = treeModel.getNode(i);
            nodeHeights[node.getNumber()] = treeModel.getNodeHeight(node);
        }
        double mappedLogLikelihood;
        try {
            mappedLogLikelihood = likelihoodDelegate.prepareStochasticMap(preparedInput, thetaBuffer(),
                    branchRateBufferOrNull(), nodeHeights, checkProbabilities, stochasticMapSampler);
        } catch (MascotLikelihoodDelegate.NumericalException e) {
            throw new IllegalStateException("MASCOT stochastic map cannot be drawn for the current parameter " +
                    "values: " + e.getMessage(), e);
        }
        if (!warnedStochasticMapLogLikelihood &&
                Math.abs(mappedLogLikelihood - logLikelihood) > 1.0e-9 * Math.max(1.0, Math.abs(logLikelihood))) {
            Logger.getLogger("dr.evomodel").warning("MASCOT stochastic map was conditioned on a forward pass with " +
                    "log likelihood " + mappedLogLikelihood + ", but the cached value is " + logLikelihood);
            warnedStochasticMapLogLikelihood = true;
        }
        stochasticMapPrepared = true;
    }

    public Parameter getMigrationRates() {
        return dynamics.getMigrationRates();
    }

    public String getMigrationGradientCompatibilityError() {
        return dynamics.getMigrationGradientCompatibilityError();
    }

    public Parameter getPopSizes() {
        return dynamics.getPopSizes();
    }

    public MascotDynamics getDynamics() {
        return dynamics;
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public double getMaxStep() {
        return maxStep;
    }

    @Override
    protected void makeDirtyInternal() {
        preparedInputKnown = false;
        likelihoodDelegateKnown = false;
    }

    @Override
    protected void invalidateDerivedState() {
        gradientKnown = false;
        ancestralStatesKnown = false;
        stochasticMapPrepared = false;
        stochasticMap = null;
    }

    @Override
    protected void handleStructuredModelChangedEvent(Model model, Object object, int index) {
        if (model == treeIntervals) {
            preparedInputKnown = false;
        } else if (!isMigrationModel(model) &&
                model != dynamics.getPopulationSizeModel() &&
                model != branchRateModel) {
            likelihoodDelegateKnown = false;
        }
    }

    @Override
    protected void handleStructuredVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == dynamics.getEpochTimes()) {
            likelihoodDelegateKnown = false;
        } else if (variable != dynamics.getMigrationRates()) {
            preparedInputKnown = false;
            likelihoodDelegateKnown = false;
        }
    }

    @Override
    protected void storeStructuredState() {
        storedDerivedIndex = currentDerivedIndex;
        storedPreparedInput = preparedInput;
        storedLikelihoodDelegate = likelihoodDelegate;
        storedPreparedInputKnown = preparedInputKnown;
        storedLikelihoodDelegateKnown = likelihoodDelegateKnown;
        storedGradientKnown = gradientKnown;
        storedAncestralStatesKnown = ancestralStatesKnown;
        stochasticMap = null;
    }

    @Override
    protected void restoreStructuredState() {
        currentDerivedIndex = storedDerivedIndex;
        preparedInput = storedPreparedInput;
        likelihoodDelegate = storedLikelihoodDelegate;
        preparedInputKnown = storedPreparedInputKnown;
        likelihoodDelegateKnown = storedLikelihoodDelegateKnown;
        gradientKnown = storedGradientKnown;
        ancestralStatesKnown = storedAncestralStatesKnown;
        stochasticMapPrepared = false;
        stochasticMap = null;
    }

    private DerivedBuffers currentDerived() {
        return derivedBuffers[currentDerivedIndex];
    }

    private void ensurePreparedInput() {
        if (tipDataCache == null) {
            tipDataCache = buildStructuredTipData(true, "tip-state attributePatterns");
        }
        if (!preparedInputKnown) {
            StructuredCoalescentSchedule schedule = StructuredCoalescentSchedule.fromTreeIntervals(
                    treeModel, treeIntervals, true, false);
            preparedInput = MascotPreparedInput.prepare(schedule, tipDataCache);
            preparedInputKnown = true;
        }
    }

    private void ensureLikelihoodDelegate() {
        if (!likelihoodDelegateKnown) {
            likelihoodDelegate = new GenericMascotLikelihoodDelegate(dynamics.getStateCount(), dynamics.getBoundaries(), maxStep);
            likelihoodDelegateKnown = true;
        }
    }

    private double[] thetaBuffer() {
        if (thetaBuffer == null) {
            thetaBuffer = new double[dynamics.getParameterCount()];
        }
        dynamics.writeThetaValues(thetaBuffer);
        return thetaBuffer;
    }

    private double[] branchRateBufferOrNull() {
        if (branchRateModel == null) {
            return null;
        }
        if (branchRateBuffer == null) {
            branchRateBuffer = new double[treeModel.getNodeCount()];
        }
        writeBranchRates(branchRateBuffer);
        return branchRateBuffer;
    }

    private void writeBranchRates(double[] destination) {
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            destination[node.getNumber()] =
                    treeModel.isRoot(node) ? 0.0 : branchRateModel.getBranchRate(treeModel, node);
        }
    }

    private void ensureDerivedBuffers(DerivedBuffers buffers, boolean needsClockGradient) {
        if (buffers.combinedGradient == null) {
            buffers.combinedGradient = new double[dynamics.getParameterCount()];
        }
        if (needsClockGradient && buffers.clockGradient == null) {
            buffers.clockGradient = new double[treeModel.getNodeCount()];
        }
    }

    private boolean isMigrationModel(Model model) {
        SubstitutionModel[] migrationModels = dynamics.getMigrationModels();
        if (migrationModels == null) {
            return false;
        }
        for (SubstitutionModel migrationModel : migrationModels) {
            if (model == migrationModel) {
                return true;
            }
        }
        return false;
    }

    @Override
    protected double calculateLogLikelihood() {
        ensurePreparedInput();
        ensureLikelihoodDelegate();
        try {
            double value = likelihoodDelegate.calculateLikelihood(preparedInput, thetaBuffer(),
                    branchRateBufferOrNull(), checkProbabilities);
            return Double.isFinite(value) ? value : Double.NEGATIVE_INFINITY;
        } catch (MascotLikelihoodDelegate.NumericalException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private void evaluateLikelihoodAndGradient(boolean failOnGradientFailure) {
        ensurePreparedInput();
        ensureLikelihoodDelegate();
        double[] branchRates = branchRateBufferOrNull();
        currentDerivedIndex = 1 - storedDerivedIndex;
        DerivedBuffers buffers = currentDerived();
        ensureDerivedBuffers(buffers, branchRates != null);

        double logLikelihood;
        try {
            logLikelihood = likelihoodDelegate.calculateLikelihoodAndDerivatives(preparedInput, thetaBuffer(), branchRates,
                    buffers.combinedGradient, buffers.clockGradient, null, checkProbabilities);
        } catch (MascotLikelihoodDelegate.NumericalException e) {
            if (failOnGradientFailure) {
                throw new IllegalStateException("MASCOT gradient cannot be evaluated for the current " +
                        "parameter values: " + e.getMessage(), e);
            }
            cacheLogLikelihood(Double.NEGATIVE_INFINITY);
            gradientKnown = false;
            return;
        }
        if (!Double.isFinite(logLikelihood)) {
            if (failOnGradientFailure) {
                throw new IllegalStateException("MASCOT gradient cannot be evaluated because the current " +
                        "log likelihood is " + logLikelihood);
            }
            cacheLogLikelihood(Double.NEGATIVE_INFINITY);
            gradientKnown = false;
            return;
        }
        if (failOnGradientFailure) {
            validateGradient(buffers.combinedGradient);
        } else if (!isValidGradient(buffers.combinedGradient)) {
            cacheLogLikelihood(logLikelihood);
            gradientKnown = false;
            return;
        }
        cacheLogLikelihood(logLikelihood);
        gradientKnown = true;
    }

    private void evaluateAncestralStates() {
        ensurePreparedInput();
        ensureLikelihoodDelegate();
        double[] branchRates = branchRateBufferOrNull();
        currentDerivedIndex = 1 - storedDerivedIndex;
        DerivedBuffers buffers = currentDerived();
        ensureDerivedBuffers(buffers, branchRates != null);
        if (buffers.ancestralStateScores == null) {
            buffers.ancestralStateScores = new double[treeModel.getNodeCount() * getStateCount()];
        }

        double logLikelihood;
        try {
            logLikelihood = likelihoodDelegate.calculateLikelihoodAndDerivatives(preparedInput, thetaBuffer(), branchRates,
                    buffers.combinedGradient, buffers.clockGradient, buffers.ancestralStateScores, checkProbabilities);
        } catch (MascotLikelihoodDelegate.NumericalException e) {
            throw new IllegalStateException("MASCOT ancestral states cannot be evaluated for the current " +
                    "parameter values: " + e.getMessage(), e);
        }
        if (!Double.isFinite(logLikelihood)) {
            throw new IllegalStateException("MASCOT ancestral states cannot be evaluated because the current " +
                    "log likelihood is " + logLikelihood);
        }
        validateGradient(buffers.combinedGradient);
        validateAncestralStateScores(buffers.ancestralStateScores);
        cacheLogLikelihood(logLikelihood);
        gradientKnown = true;
        ancestralStatesKnown = true;
    }

    /**
     * Checks that every internal-node row was written and is finite. Row sums are not checked: row entry s is
     * P[s](A[s] - A.P + a_l) with a_l = 1 and sum(P) = 1, so the row sums to 1 for any adjoint A. A sum test can
     * detect only roundoff, never discretization or adjoint errors.
     */
    private void validateAncestralStateScores(double[] scores) {
        if (scores == null) {
            throw new IllegalStateException("MASCOT ancestral-state evaluation returned no scores");
        }
        int stateCount = getStateCount();
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (treeModel.isExternal(node)) {
                continue;
            }
            int nodeNumber = node.getNumber();
            int offset = nodeNumber * stateCount;
            for (int s = 0; s < stateCount; s++) {
                double v = scores[offset + s];
                if (Double.isNaN(v)) {
                    throw new IllegalStateException("MASCOT ancestral-state row was never written for " +
                            "internal node " + nodeNumber + " (state " + s + " is NaN)");
                }
                if (!Double.isFinite(v)) {
                    throw new IllegalStateException("MASCOT ancestral-state score is non-finite at node " +
                            nodeNumber + ", state " + s + ": " + v);
                }
            }
        }
    }

    private static void validateGradient(double[] gradient) {
        if (gradient == null) {
            throw new IllegalStateException("MASCOT gradient evaluation returned no gradient");
        }
        for (int i = 0; i < gradient.length; i++) {
            if (!Double.isFinite(gradient[i])) {
                throw new IllegalStateException("MASCOT gradient contains a non-finite entry at index " +
                        i + ": " + gradient[i]);
            }
        }
    }

    private static boolean isValidGradient(double[] gradient) {
        if (gradient == null) {
            return false;
        }
        for (int i = 0; i < gradient.length; i++) {
            if (!Double.isFinite(gradient[i])) {
                return false;
            }
        }
        return true;
    }

    public double[] getGradientLogDensity(StructuredCoalescentLikelihoodGradient wrt) {
        switch (wrt.getType()) {
            case MIGRATION_RATE:
                return getMigrationGradientLogDensity();
            case POPULATION_SIZE:
                return getPopSizeGradientLogDensity();
            case CLOCK_RATE:
                return getClockRateGradientLogDensity();
            default:
                throw new IllegalArgumentException("Unsupported wrtParameter for a MASCOT likelihood: " + wrt.getType());
        }
    }

    public double[] getClockRateGradientLogDensity() {
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new IllegalStateException("clockRate gradient requires a branchRateModel that implements " +
                    "DifferentiableBranchRates, got: " +
                    (branchRateModel == null ? "null (no branchRateModel supplied)" : branchRateModel.getClass().getName()));
        }
        DifferentiableBranchRates branchRates = (DifferentiableBranchRates) branchRateModel;
        double[] raw = getClockGradientLogDensity();
        double[] result = new double[treeModel.getNodeCount() - 1];
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
                result[branchRates.getParameterIndexFromNode(node)] = raw[node.getNumber()];
            }
        }
        return branchRates.updateGradientLogDensity(result, null, 0, result.length);
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return "";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.emptyList();
    }
}
