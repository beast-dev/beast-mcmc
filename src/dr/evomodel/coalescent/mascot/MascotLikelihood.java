/*
 * MascotLikelihood.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.PatternList;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

/**
 * BEAST-X model wrapper for {@link MascotCore}.
 */
public class MascotLikelihood extends AbstractModelLikelihood implements Reportable {

    public static final String MASCOT_LIKELIHOOD = "mascotLikelihood";

    private final TreeModel treeModel;
    // Shared with BASTA's BastaLikelihood: BigFastTreeIntervals is the single
    // source of truth for "walk this tree's coalescent/sample events in time
    // order," instead of MascotEventTape re-deriving that order itself.
    private final BestSignalsFromBigFastTreeIntervals treeIntervals;
    private final PatternList tipPatterns;
    private final MascotDynamics dynamics;
    private final double maxStep;
    private final boolean checkProbabilities;
    // Nullable: null means no clock scaling (branchRates=null passed into
    // MascotCore.evaluate(...), the exact same code path as before this field
    // existed -- see MascotCore's own null-branchRates documentation).
    private final BranchRateModel branchRateModel;

    private MascotEventTape eventTape;
    private MascotCore core;
    private boolean eventTapeKnown;
    private boolean coreKnown;
    private boolean likelihoodKnown;
    private boolean gradientKnown;
    private boolean gradientPrimedLikelihood;

    private double logLikelihood;
    // The combined, MascotCore-native flat gradient (epoch-major [migration, popSizes]).
    // MascotGradient slices this per-part; callers outside this package never see it directly.
    private double[] combinedGradient;
    // d(logLikelihood)/d(branchRate[lineageId]); null unless branchRateModel != null.
    private double[] clockGradient;

    // Snapshots taken by storeState() and restored by restoreState() on a
    // rejected proposal. Restoring the actual eventTape/core references (rather
    // than just invalidating eventTapeKnown/coreKnown, as an earlier version of
    // this class did) means a rejected proposal that never touched the tree or
    // epoch times does not force a full event-tape/core rebuild on the next
    // evaluation. Both MascotEventTape and MascotCore are effectively immutable
    // from this wrapper's point of view once built (MascotCore's internal
    // per-epoch rate cache is fully recomputed from theta at the start of every
    // evaluate() call regardless), so sharing the stored reference back in is safe.
    private MascotEventTape storedEventTape;
    private MascotCore storedCore;
    private boolean storedEventTapeKnown;
    private boolean storedCoreKnown;
    private boolean storedLikelihoodKnown;
    private boolean storedGradientKnown;
    private double storedLogLikelihood;
    private double[] storedCombinedGradient;
    private double[] storedClockGradient;

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            Parameter migrationRates,
                            Parameter popSizes,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities) {
        this(name, treeModel, tipPatterns, migrationRates, popSizes, epochTimes, stateCount, maxStep,
                checkProbabilities, null);
    }

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            Parameter migrationRates,
                            Parameter popSizes,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities,
                            BranchRateModel branchRateModel) {
        this(name, treeModel, tipPatterns, migrationRates, null, popSizes, epochTimes, stateCount, maxStep,
                checkProbabilities, branchRateModel);
    }

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            SubstitutionModel migrationModel,
                            Parameter popSizes,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities) {
        this(name, treeModel, tipPatterns, null, migrationModel, popSizes, epochTimes, stateCount, maxStep,
                checkProbabilities, null);
    }

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            SubstitutionModel migrationModel,
                            Parameter popSizes,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities,
                            BranchRateModel branchRateModel) {
        this(name, treeModel, tipPatterns, null, migrationModel, popSizes, epochTimes, stateCount, maxStep,
                checkProbabilities, branchRateModel);
    }

    private MascotLikelihood(String name,
                            TreeModel treeModel,
                            PatternList tipPatterns,
                            Parameter migrationRates,
                            SubstitutionModel migrationModel,
                            Parameter popSizes,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities,
                            BranchRateModel branchRateModel) {
        super(name == null ? MASCOT_LIKELIHOOD : name);
        this.treeModel = treeModel;
        this.treeIntervals = new BestSignalsFromBigFastTreeIntervals(treeModel);
        this.tipPatterns = tipPatterns;
        this.dynamics = migrationRates != null ?
                new MascotDynamics(stateCount, migrationRates, popSizes, epochTimes) :
                new MascotDynamics(stateCount, migrationModel, popSizes, epochTimes);
        this.maxStep = maxStep;
        this.checkProbabilities = checkProbabilities;
        this.branchRateModel = branchRateModel;

        addModel(treeIntervals);
        // tipPatterns is fixed input data (like a sequence alignment), not an
        // estimated model variable, so it is not registered as a listened
        // variable here -- it never changes over the course of an analysis.
        if (migrationRates != null) {
            addVariable(migrationRates);
        } else {
            addModel(migrationModel);
        }
        addVariable(popSizes);
        if (epochTimes != null) {
            addVariable(epochTimes);
        }
        if (branchRateModel != null) {
            addModel(branchRateModel);
        }

        this.eventTapeKnown = false;
        this.coreKnown = false;
        this.likelihoodKnown = false;
        this.gradientKnown = false;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            if (gradientPrimedLikelihood && !gradientKnown) {
                evaluateLikelihoodAndGradient(false);
            } else {
                evaluateLikelihoodOnly();
            }
        }
        return logLikelihood;
    }

    /**
     * MascotCore's own combined flat gradient (epoch-major [migration, popSizes]).
     * Most callers want {@link #getMigrationGradientLogDensity()} or
     * {@link #getPopSizeGradientLogDensity()} instead (see {@link MascotGradient}).
     */
    public double[] getGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        // Clone on every return, not just on a cache hit: the stored gradient array
        // is reused across calls, so a caller that mutated the returned array in
        // place would otherwise silently corrupt the cache for the next caller.
        return combinedGradient.clone();
    }

    public double[] getMigrationGradientLogDensity() {
        return dynamics.extractMigrationGradient(getGradientLogDensity());
    }

    public double[] getPopSizeGradientLogDensity() {
        return dynamics.extractPopSizeGradient(getGradientLogDensity());
    }

    /**
     * d(logLikelihood)/d(branchRate[lineageId]), indexed by tree node number (see
     * {@link MascotEventTape#buildBranchRates}). Only valid when {@link
     * #getBranchRateModel()} is non-null.
     */
    public double[] getClockGradientLogDensity() {
        if (!gradientKnown) {
            evaluateLikelihoodAndGradient(true);
        }
        return clockGradient == null ? null : clockGradient.clone();
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

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public int getStateCount() {
        return dynamics.getStateCount();
    }

    public double getMaxStep() {
        return maxStep;
    }

    void enableGradientPrimedLikelihood() {
        gradientPrimedLikelihood = true;
    }

    @Override
    public void makeDirty() {
        eventTapeKnown = false;
        coreKnown = false;
        likelihoodKnown = false;
        gradientKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeIntervals) {
            eventTapeKnown = false;
        } else if (model == dynamics.getMigrationModel()) {
            // A substitution-model-backed migration-rate change only changes the
            // numeric Q matrix read into theta; event tape and epoch boundaries
            // remain valid.
        } else if (model == branchRateModel) {
            // A clock-rate parameter changed. The event tape and MascotCore's
            // per-epoch rate cache are both keyed on the tree/theta, neither of
            // which changed here -- only the numeric likelihood/gradient are
            // stale. Falling through to the generic "else" below would set
            // coreKnown=false and force a full MascotCore/tape rebuild on every
            // HMC step touching the clock rate, defeating the tape-reuse design
            // (see MascotCore's OperationTapeStore doc comment).
        } else {
            coreKnown = false;
        }
        likelihoodKnown = false;
        gradientKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == dynamics.getMigrationRates() || variable == dynamics.getPopSizes()) {
            // The event tape is still valid; only the numeric likelihood and gradient change.
        } else if (variable == dynamics.getEpochTimes()) {
            coreKnown = false;
        } else {
            eventTapeKnown = false;
            coreKnown = false;
        }
        likelihoodKnown = false;
        gradientKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedCombinedGradient = combinedGradient;
        storedClockGradient = clockGradient;
        storedEventTape = eventTape;
        storedCore = core;
        storedEventTapeKnown = eventTapeKnown;
        storedCoreKnown = coreKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedGradientKnown = gradientKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        combinedGradient = storedCombinedGradient;
        clockGradient = storedClockGradient;
        eventTape = storedEventTape;
        core = storedCore;
        eventTapeKnown = storedEventTapeKnown;
        coreKnown = storedCoreKnown;
        likelihoodKnown = storedLikelihoodKnown;
        gradientKnown = storedGradientKnown;
    }

    @Override
    protected void acceptState() {
        // Nothing to do: current fields already reflect the accepted state.
    }

    @Override
    public String getReport() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    private void ensureEventTape() {
        if (!eventTapeKnown) {
            eventTape = MascotEventTape.fromTree(treeIntervals, tipPatterns, dynamics.getStateCount());
            eventTapeKnown = true;
        }
    }

    private void ensureCore() {
        if (!coreKnown) {
            core = new MascotCore(dynamics.getStateCount(), dynamics.getBoundaries(), maxStep);
            coreKnown = true;
        }
    }

    private double[] branchRatesOrNull() {
        return branchRateModel == null ? null : MascotEventTape.buildBranchRates(treeModel, branchRateModel);
    }

    private void evaluateLikelihoodOnly() {
        ensureEventTape();
        ensureCore();
        try {
            logLikelihood = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(),
                    branchRatesOrNull(), false, checkProbabilities, false).logLikelihood;
            if (!Double.isFinite(logLikelihood)) {
                logLikelihood = Double.NEGATIVE_INFINITY;
            }
        } catch (MascotCore.NumericalException e) {
            logLikelihood = Double.NEGATIVE_INFINITY;
        }
        likelihoodKnown = true;
    }

    private void evaluateLikelihoodAndGradient(boolean failOnGradientFailure) {
        ensureEventTape();
        ensureCore();
        MascotCore.Result result;
        try {
            result = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(),
                    branchRatesOrNull(), true, checkProbabilities, false);
        } catch (MascotCore.NumericalException e) {
            if (failOnGradientFailure) {
                throw new IllegalStateException("MASCOT gradient cannot be evaluated for the current " +
                        "parameter values: " + e.getMessage(), e);
            }
            logLikelihood = Double.NEGATIVE_INFINITY;
            likelihoodKnown = true;
            gradientKnown = false;
            return;
        }
        if (!Double.isFinite(result.logLikelihood)) {
            if (failOnGradientFailure) {
                throw new IllegalStateException("MASCOT gradient cannot be evaluated because the current " +
                        "log likelihood is " + result.logLikelihood);
            }
            logLikelihood = Double.NEGATIVE_INFINITY;
            likelihoodKnown = true;
            gradientKnown = false;
            return;
        }
        if (failOnGradientFailure) {
            validateGradient(result.gradient);
        } else if (!isValidGradient(result.gradient)) {
            logLikelihood = result.logLikelihood;
            likelihoodKnown = true;
            gradientKnown = false;
            return;
        }
        logLikelihood = result.logLikelihood;
        combinedGradient = result.gradient;
        clockGradient = result.clockGradient;
        likelihoodKnown = true;
        gradientKnown = true;
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
}
