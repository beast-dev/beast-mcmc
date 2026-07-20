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
import dr.evomodel.coalescent.AbstractStructuredCoalescentLikelihood;
import dr.evomodel.coalescent.StructuredTipStates;
import dr.evomodel.coalescent.basta.AbstractPopulationSizeModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * @author Filippo Monti
 */
public class MascotLikelihood extends AbstractStructuredCoalescentLikelihood
        implements Citable {

    public static final String MASCOT_LIKELIHOOD = "mascotLikelihood";

    /**
     * Default {@link TreeTrait} name under which adjoint node-state scores
     * (see {@link #getAncestralStateScores()}) are exposed to tree logging,
     * used when no explicit tag name is given to the constructor (see
     * {@code stateTagName} in {@code StructuredCoalescentLikelihoodParser}).
     * Named {@code stateSensitivity}, not {@code states}: these are exact
     * sensitivities of the discretized RK4 likelihood
     * (d(logL)/d(alpha_v,s)), not a posterior ancestral-state
     * reconstruction -- a nonnegativity survey found a converged,
     * step-size-independent negative value under highly unequal population
     * sizes, so this quantity is not always a valid probability
     * distribution (see MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md's
     * "Probability interpretation" note) and must not be presented or named
     * as one. Also deliberately distinct from BASTA's "states" tag, which is
     * a genuine up-down MAP/joint reconstruction -- the two should not look
     * interchangeable in output by default.
     */
    public static final String DEFAULT_ANCESTRAL_STATE_TAG_NAME = "mascot.stateSensitivity";

    private final TreeModel treeModel;
    private final PatternList tipPatterns;
    private final MascotDynamics dynamics;
    private final double maxStep;
    private final boolean checkProbabilities;
    private final String ancestralStateTagName;

    private MascotEventTape eventTape;
    private MascotCore core;
    private boolean eventTapeKnown;
    private boolean coreKnown;
    private boolean gradientKnown;
    // Tip states are fixed input data, never re-estimated, so unlike eventTape
    // (rebuilt on every tree change) this is built once and never invalidated.
    // Avoids re-deriving every tip's partial vector -- an O(tipCount) taxon-name
    // lookup plus a fresh allocation, per tip -- on every tree-changing operator.
    private double[][] tipPartialsCache;

    // The combined, MascotCore-native flat gradient (epoch-major [migration, popSizes]).
    // MascotGradient slices this per-part; callers outside this package never see it directly.
    private double[] combinedGradient;
    // d(logLikelihood)/d(branchRate[lineageId]); null unless branchRateModel != null.
    private double[] clockGradient;
    private boolean ancestralStatesKnown;
    // Flat, node-major adjoint node-state scores (see MascotCore.Result#ancestralStateScores);
    // null until getAncestralStateScores() is first called.
    private double[] ancestralStateScores;

    // Snapshots taken by storeStructuredState() and restored by restoreStructuredState()
    // on a rejected proposal. Restoring the actual eventTape/core references (rather
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
    private boolean storedGradientKnown;
    private boolean storedAncestralStatesKnown;
    private double[] storedCombinedGradient;
    private double[] storedClockGradient;
    private double[] storedAncestralStateScores;

    /**
     * Exactly one of {@code migrationRates} (a raw parameter of native
     * positive migration rates) or {@code migrationModels} (one
     * {@link SubstitutionModel} per epoch -- a single-element array is the
     * constant-through-time case; see {@link MascotDynamics}'s class doc) must
     * be non-null; the other must be null. {@code branchRateModel} may be
     * null (no clock gradient support). Callers that don't need a custom
     * ancestral-state tree-trait name should pass {@link
     * #DEFAULT_ANCESTRAL_STATE_TAG_NAME}. This is the one constructor:
     * argument defaulting (e.g. {@code null} branchRateModel, {@code
     * DEFAULT_ANCESTRAL_STATE_TAG_NAME}) is the XML parser's job (see {@code
     * StructuredCoalescentLikelihoodParser.parseMascot}), not a further stack
     * of overloads here.
     */
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
                            String ancestralStateTagName) {
        super(name == null ? MASCOT_LIKELIHOOD : name, treeModel, tipPatterns, stateCount, branchRateModel);
        this.treeModel = treeModel;
        this.tipPatterns = tipPatterns;
        this.dynamics = migrationRates != null ?
                new MascotDynamics(stateCount, migrationRates, populationSizeModel, epochTimes) :
                new MascotDynamics(stateCount, migrationModels, populationSizeModel, epochTimes);
        this.maxStep = maxStep;
        this.checkProbabilities = checkProbabilities;
        this.ancestralStateTagName = ancestralStateTagName;

        // tipPatterns is fixed input data (like a sequence alignment), not an
        // estimated model variable, so it is not registered as a listened
        // variable here -- it never changes over the course of an analysis.
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

        this.eventTapeKnown = false;
        this.coreKnown = false;
        this.gradientKnown = false;

        // TreeTraitProvider.Helper keys a trait by its getTraitName() at
        // registration time (not read live afterward), so the tag name is
        // captured once here from the immutable ancestralStateTagName
        // constructor parameter (like BASTA's own immutable "tag" field),
        // not a post-construction setter.
        treeTraits.addTrait(new TreeTrait.DA() {
            public String getTraitName() {
                return ancestralStateTagName;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            // getTraitClass() is inherited from TreeTrait.DA (double[].class);
            // do not override it -- an earlier version incorrectly overrode it
            // with Double.class.

            public double[] getTrait(Tree tree, NodeRef node) {
                // Tip rows are intentionally left NaN by MascotCore (Section
                // 5.3 of the design doc); generic tree-trait logging (Intent.NODE,
                // the default ALL restriction) still requests every node's
                // value, tips included, so return the observed/uncertain tip
                // partial here rather than an all-NaN vector.
                if (tree.isExternal(node)) {
                    return StructuredTipStates.getPartials(tree, node, tipPatterns, getStateCount(), true,
                            ancestralStateTagName + " tip trait");
                }
                double[] row = new double[getStateCount()];
                getAncestralStateScores(node.getNumber(), row);
                return row;
            }
        });
    }

    // getLogLikelihood() itself, getTreeTraits()/getTreeTrait(), and
    // getBranchRateModel() are inherited unchanged from
    // AbstractStructuredCoalescentLikelihood: calculateLogLikelihood() below
    // is the only hook it calls, and it deliberately never opportunistically
    // computes the gradient there, even if some HMC operator elsewhere primed
    // this likelihood -- an operator that only needs the likelihood (tree
    // moves, bitFlip, plain scaleOperator) must not pay for the adjoint pass
    // on every one of its evaluations just because an unrelated
    // structuredCoalescentLikelihoodGradient element also exists in the XML.
    // getGradientLogDensity() below is the only path that ever computes the
    // gradient, and only when actually called; when it is, it caches the
    // likelihood too (via cacheLogLikelihood()), so a getLogLikelihood() call
    // right after a real HMC step is still a cache hit.

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

    /**
     * Flat, node-major adjoint node-state scores for every internal tree
     * node, indexed by {@code nodeNumber * getStateCount() + state} (see
     * {@code MascotCore.Result#ancestralStateScores} and
     * {@code MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md}). Tip rows and any
     * unused node-number slots are {@link Double#NaN}. These are exact
     * sensitivities of the discretized RK4 likelihood -- d(logL)/d(alpha_v,s),
     * not a posterior-probability interpretation: a broad nonnegativity
     * survey found a converged, step-size-independent negative score under
     * highly unequal population sizes, so this quantity is not always a
     * valid probability distribution (see the design document's
     * "Probability interpretation" note under Section 19), and legitimate
     * negative entries are permitted here, not treated as failures. This
     * call still throws on a genuine correctness violation: a non-finite
     * entry, a missing (never-written) internal-node row, or a row that
     * doesn't sum to one within tolerance (an identity that holds regardless
     * of individual entries' sign, so its violation always indicates a bug,
     * never a valid parameter regime).
     */
    public double[] getAncestralStateScores() {
        if (!ancestralStatesKnown) {
            evaluateAncestralStates();
        }
        return ancestralStateScores.clone();
    }

    /** Row-copy convenience for allocation-conscious tree logging; avoids exposing the flat layout to callers. */
    public void getAncestralStateScores(int nodeNumber, double[] destination) {
        if (destination.length != getStateCount()) {
            throw new IllegalArgumentException("destination length " + destination.length +
                    " does not match stateCount " + getStateCount());
        }
        if (!ancestralStatesKnown) {
            evaluateAncestralStates();
        }
        if (nodeNumber < 0 || (nodeNumber + 1) * getStateCount() > ancestralStateScores.length) {
            throw new IllegalArgumentException("node number out of range: " + nodeNumber);
        }
        System.arraycopy(ancestralStateScores, nodeNumber * getStateCount(), destination, 0, getStateCount());
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
        eventTapeKnown = false;
        coreKnown = false;
    }

    @Override
    protected void invalidateDerivedState() {
        gradientKnown = false;
        ancestralStatesKnown = false;
    }

    // Base class calls this, then unconditionally invalidates the shared
    // likelihood cache and derived state (invalidateDerivedState() above),
    // then fires the model-changed event -- so this only needs the
    // classification logic for eventTapeKnown/coreKnown.
    @Override
    protected void handleStructuredModelChangedEvent(Model model, Object object, int index) {
        if (model == treeIntervals) {
            eventTapeKnown = false;
        } else if (isMigrationModel(model)) {
            // A substitution-model-backed migration-rate change (in any
            // epoch) only changes the numeric Q matrix read into theta;
            // event tape and epoch boundaries remain valid.
        } else if (model == dynamics.getPopulationSizeModel()) {
            // A population-size-model value change (Constant or
            // PiecewiseConstant) only changes the numeric logN values read
            // into theta; event tape and epoch boundaries remain valid.
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
    }

    @Override
    protected void handleStructuredVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == dynamics.getMigrationRates()) {
            // The event tape is still valid; only the numeric likelihood and gradient change.
        } else if (variable == dynamics.getEpochTimes()) {
            coreKnown = false;
        } else {
            eventTapeKnown = false;
            coreKnown = false;
        }
    }

    @Override
    protected void storeStructuredState() {
        storedCombinedGradient = combinedGradient;
        storedClockGradient = clockGradient;
        storedAncestralStateScores = ancestralStateScores;
        storedEventTape = eventTape;
        storedCore = core;
        storedEventTapeKnown = eventTapeKnown;
        storedCoreKnown = coreKnown;
        storedGradientKnown = gradientKnown;
        storedAncestralStatesKnown = ancestralStatesKnown;
    }

    @Override
    protected void restoreStructuredState() {
        combinedGradient = storedCombinedGradient;
        clockGradient = storedClockGradient;
        ancestralStateScores = storedAncestralStateScores;
        eventTape = storedEventTape;
        core = storedCore;
        eventTapeKnown = storedEventTapeKnown;
        coreKnown = storedCoreKnown;
        gradientKnown = storedGradientKnown;
        ancestralStatesKnown = storedAncestralStatesKnown;
    }

    // acceptState() and getReport() are inherited unchanged: both engines'
    // acceptState() was already a no-op, and this class's getReport() always
    // matched AbstractStructuredCoalescentLikelihood.getDefaultReport() exactly.

    private void ensureEventTape() {
        // validateSinglePattern(tipPatterns, stateCount, ...) already ran in
        // AbstractStructuredCoalescentLikelihood's constructor, shared with BASTA.
        if (tipPartialsCache == null) {
            tipPartialsCache = StructuredTipStates.buildPartialsCache(treeModel, tipPatterns, dynamics.getStateCount(),
                    true, "tip-state attributePatterns");
        }
        if (!eventTapeKnown) {
            eventTape = MascotEventTape.fromTree(treeIntervals, tipPartialsCache, dynamics.getStateCount());
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
        ensureEventTape();
        ensureCore();
        try {
            double value = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(),
                    branchRatesOrNull(), false, checkProbabilities, false).logLikelihood;
            return Double.isFinite(value) ? value : Double.NEGATIVE_INFINITY;
        } catch (MascotCore.NumericalException e) {
            return Double.NEGATIVE_INFINITY;
        }
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
            cacheLogLikelihood(Double.NEGATIVE_INFINITY);
            gradientKnown = false;
            return;
        }
        if (!Double.isFinite(result.logLikelihood)) {
            if (failOnGradientFailure) {
                throw new IllegalStateException("MASCOT gradient cannot be evaluated because the current " +
                        "log likelihood is " + result.logLikelihood);
            }
            cacheLogLikelihood(Double.NEGATIVE_INFINITY);
            gradientKnown = false;
            return;
        }
        if (failOnGradientFailure) {
            validateGradient(result.gradient);
        } else if (!isValidGradient(result.gradient)) {
            cacheLogLikelihood(result.logLikelihood);
            gradientKnown = false;
            return;
        }
        cacheLogLikelihood(result.logLikelihood);
        combinedGradient = result.gradient;
        clockGradient = result.clockGradient;
        gradientKnown = true;
    }

    /**
     * Evaluates the likelihood together with adjoint node-state scores for
     * every internal tree node. Also requests the parameter gradient from
     * the same MascotCore evaluation (rather than a separate one) so it is
     * cached "for free" alongside the ancestral scores when both are needed;
     * this does not run more often than plain likelihood/gradient evaluation
     * since it is only invoked by {@link #getAncestralStateScores()} and its
     * row-copy overload, which tree logging calls far less frequently than
     * every MCMC/HMC step.
     */
    private void evaluateAncestralStates() {
        ensureEventTape();
        ensureCore();
        MascotCore.Result result;
        try {
            result = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(),
                    branchRatesOrNull(), true, true, checkProbabilities, false);
        } catch (MascotCore.NumericalException e) {
            throw new IllegalStateException("MASCOT ancestral states cannot be evaluated for the current " +
                    "parameter values: " + e.getMessage(), e);
        }
        if (!Double.isFinite(result.logLikelihood)) {
            throw new IllegalStateException("MASCOT ancestral states cannot be evaluated because the current " +
                    "log likelihood is " + result.logLikelihood);
        }
        validateGradient(result.gradient);
        validateAncestralStateScores(result.ancestralStateScores);
        cacheLogLikelihood(result.logLikelihood);
        combinedGradient = result.gradient;
        clockGradient = result.clockGradient;
        ancestralStateScores = result.ancestralStateScores;
        gradientKnown = true;
        ancestralStatesKnown = true;
    }

    /**
     * Checks finiteness, that every internal node's row was actually
     * written, and sum-to-one for every internal-node row. Iterates the
     * real tree (rather than inferring tip-vs-internal from whether a row's
     * first entry happens to be NaN) so a genuinely missing internal-node
     * row -- e.g. from a tape/indexing bug -- is caught as a failure
     * instead of silently treated like a tip.
     * <p/>
     * Deliberately does <em>not</em> reject negative entries: this quantity
     * is a sensitivity, not a probability (see {@link #getAncestralStateScores()}),
     * and a broad nonnegativity survey already found legitimate, converged
     * negative values under highly unequal population sizes -- rejecting
     * them here would mean a valid MCMC/HMC state could make tree logging
     * throw. The sum-to-one check is still a hard failure: that identity
     * (Section 3 of the design doc) holds regardless of individual entries'
     * sign, so its violation always indicates a bug, never a valid regime.
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
            double sum = 0.0;
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
                sum += v;
            }
            if (Math.abs(sum - 1.0) > 1.0e-8) {
                throw new IllegalStateException("MASCOT ancestral-state row sum " + sum +
                        " does not equal 1 within tolerance 1e-8 at node " + nodeNumber + " (maxStep=" + maxStep + ")");
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
        return List.of();
    }
}
