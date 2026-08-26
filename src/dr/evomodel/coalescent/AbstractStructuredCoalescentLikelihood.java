/*
 * AbstractStructuredCoalescentLikelihood.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.HiddenCodons;
import dr.evolution.datatype.HiddenDataType;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.basta.AbstractPopulationSizeModel;
import dr.evomodel.coalescent.basta.PiecewiseConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.StructuredCoalescentLikelihoodGradient;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.function.Function;

/**
 * @author Filippo Monti
 * @author Yucai Shao
 * @author Marc A. Suchard
 * @author Guy Baele
 *
 * Shared wrapper infrastructure for structured-coalescent likelihoods:
 * BASTA's matrix-exponential/BEAGLE engine ({@code BastaLikelihood}) and
 * MASCOT's RK4/adjoint ODE engine ({@code MascotLikelihood}). This class
 * owns the common likelihood lifecycle -- tree/pattern/state-count
 * bookkeeping, tree-trait registration, the lazy logLikelihood cache, and
 * the store/restore/dirty/model-changed plumbing every {@code
 * AbstractModelLikelihood} needs -- via template methods that call into a
 * small set of engine-specific hooks. Everything downstream of that --
 * interval traversal, transition-matrix/ODE evaluation, gradients, and
 * ancestral-state reconstruction -- is genuinely different between the two
 * engines and stays in each subclass; this class deliberately does not
 * become an abstraction for "a structured-coalescent numerical algorithm,"
 * only for "a cached likelihood over a tree, tip-state patterns, and
 * optional branch rates."
 */
public abstract class AbstractStructuredCoalescentLikelihood extends AbstractModelLikelihood
        implements TreeTraitProvider, Reportable {

    protected final Tree tree;
    protected final PatternList patternList;
    protected final int stateCount;

    protected final BestSignalsFromBigFastTreeIntervals treeIntervals;
    // Nullable: not every structured-coalescent likelihood requires a clock
    // (MASCOT allows null, meaning no rate scaling of the migration process).
    protected final BranchRateModel branchRateModel;

    protected final SubstitutionModel[] substitutionModels;

    protected final AbstractPopulationSizeModel populationSizeModel;

    protected final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper();

    // Ancestral/mode-state-label plumbing shared by BASTA's up-down
    // reconstruction and MASCOT's argmax-over-adjoint-scores mode state (see
    // addModeStateTrait below): the "index -> state name" step is identical
    // between the two engines even though how each gets its winning index is
    // not.
    protected final DataType dataType;
    protected final String tag;
    private final CodeFormatter formatter;

    protected double logLikelihood;
    protected boolean likelihoodKnown;

    private double storedLogLikelihood;
    private boolean storedLikelihoodKnown;

    protected AbstractStructuredCoalescentLikelihood(String name,
                                                     Tree tree,
                                                     PatternList patternList,
                                                     int stateCount,
                                                     BranchRateModel branchRateModel,
                                                     SubstitutionModel substitutionModel,
                                                     AbstractPopulationSizeModel populationSizeModel,
                                                     DataType dataType,
                                                     String tag) {
        this(name, tree, patternList, stateCount, branchRateModel, new SubstitutionModel[]{substitutionModel}, populationSizeModel, dataType, tag);
    }

    protected AbstractStructuredCoalescentLikelihood(String name,
                                                     Tree tree,
                                                     PatternList patternList,
                                                     int stateCount,
                                                     BranchRateModel branchRateModel,
                                                     SubstitutionModel[] substitutionModels,
                                                     AbstractPopulationSizeModel populationSizeModel,
                                                     DataType dataType,
                                                     String tag) {
        super(name);

        if (tree == null) {
            throw new IllegalArgumentException("tree must not be null");
        }
        if (patternList == null) {
            throw new IllegalArgumentException("patternList must not be null");
        }
        if (stateCount <= 0) {
            throw new IllegalArgumentException("stateCount must be positive: " + stateCount);
        }
        StructuredTipStates.validateSinglePattern(patternList, stateCount, "structured-coalescent tip-state attributePatterns");

        this.tree = tree;
        this.patternList = patternList;
        this.stateCount = stateCount;

        this.treeIntervals = buildTreeIntervals(tree);
        addModel(treeIntervals);

        this.branchRateModel = branchRateModel;
        if (branchRateModel != null) {
            addModel(branchRateModel);
        }
        this.substitutionModels = substitutionModels;

        this.populationSizeModel = populationSizeModel;
        addModel(populationSizeModel);
        if (populationSizeModel instanceof PiecewiseConstantPopulationSizeModel) {
            ((PiecewiseConstantPopulationSizeModel) populationSizeModel).setTreeIntervals(treeIntervals);
        }

        this.dataType = dataType;
        this.tag = tag;
        this.formatter = new CodeFormatter(dataType, false);

        likelihoodKnown = false;
    }

    private static BestSignalsFromBigFastTreeIntervals buildTreeIntervals(Tree tree) {
        if (tree instanceof TreeModel) {
            return new BestSignalsFromBigFastTreeIntervals((TreeModel) tree);
        }
        if (tree instanceof MutableTreeModel) {
            return new BestSignalsFromBigFastTreeIntervals("intervals", (MutableTreeModel) tree);
        }
        throw new IllegalArgumentException("Unsupported tree implementation: " + tree.getClass().getName());
    }

    @Override
    public final Model getModel() {
        return this;
    }

    public final Tree getTree() {
        return tree;
    }

    public final PatternList getPatternList() {
        return patternList;
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getPatternCount() {
        return patternList.getPatternCount();
    }

    public int getTipCount() {
        return tree.getExternalNodeCount();
    }

    protected final StructuredCoalescentTipData buildStructuredTipData(boolean allowAmbiguities, String context) {
        return StructuredCoalescentTipData.fromPartials(tree.getNodeCount(), stateCount,
                buildStructuredTipPartials(allowAmbiguities, context));
    }

    protected final double[][] buildStructuredTipPartials(boolean allowAmbiguities, String context) {
        return StructuredTipStates.buildPartialsCache(tree, patternList, stateCount, allowAmbiguities, context);
    }

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public void addTrait(TreeTrait trait) {
        treeTraits.addTrait(trait);
    }

    public void addTraits(TreeTrait[] traits) {
        treeTraits.addTraits(traits);
    }

    public final DataType getDataType() {
        return dataType;
    }

    public final String getTag() {
        return tag;
    }

    /**
     * Quotes and formats a single resolved state index through this
     * likelihood's {@link DataType}, e.g. {@code "Region1"} -- the same
     * label format {@code location="Region1"}-style tree annotation expects.
     */
    protected final String formattedState(int state) {
        formatter.reset();
        return "\"" + formatter.getCodeString(state) + "\"";
    }

    /** Index of the largest entry; ties favor the earliest index. */
    protected static int argmax(double[] values) {
        double max = values[0];
        int choice = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
                choice = i;
            }
        }
        return choice;
    }

    @FunctionalInterface
    protected interface NodeStateFunction {
        int getState(Tree tree, NodeRef node);
    }

    /**
     * Registers a {@code TreeTrait<Integer>} under {@code traitName} whose
     * {@link TreeTrait#getTraitString} is the resolved state's name via this
     * likelihood's {@link DataType} (e.g. {@code location="Region1"}), while
     * {@link TreeTrait#getTrait} stays the plain state index for programmatic
     * use. This is the plumbing BASTA's up-down reconstruction and MASCOT's
     * argmax-over-adjoint-scores mode state both need identically; only "how
     * do I get a winning state index for this node" differs per engine,
     * supplied here as {@code stateFunction}.
     */
    protected final void addModeStateTrait(String traitName, NodeStateFunction stateFunction) {
        treeTraits.addTrait(new TreeTrait.I() {
            public String getTraitName() {
                return traitName;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Integer getTrait(Tree tree, NodeRef node) {
                return stateFunction.getState(tree, node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formattedState(stateFunction.getState(tree, node));
            }
        });
    }

    /**
     * Maps a state index to its {@link DataType} code string, matching
     * {@code GeneralDataType}/{@code HiddenCodons}/{@code HiddenDataType}
     * conventions. Moved here (from what was originally a BASTA-only private
     * helper) since MASCOT's mode-state label needs exactly the same mapping.
     */
    private static final class CodeFormatter {
        private final Function<String, String> appender;
        private final Function<Integer, String> getter;
        private boolean first = true;

        CodeFormatter(DataType dataType, boolean stripHiddenState) {
            this.appender = (dataType instanceof GeneralDataType) ?
                    (codeString) -> codeString + " " : Function.identity();

            if (dataType instanceof HiddenCodons) {
                this.getter = (stripHiddenState) ?
                        ((HiddenCodons) dataType)::getTripletWithoutHiddenCode :
                        dataType::getTriplet;
            } else if (dataType instanceof HiddenDataType && stripHiddenState) {
                this.getter = ((HiddenDataType) dataType)::getCodeWithoutHiddenState;
            } else {
                this.getter = dataType::getCode;
            }
        }

        String getCodeString(int state) {
            String code = getter.apply(state);
            if (first) {
                first = false;
            } else {
                code = appender.apply(code);
            }
            return code;
        }

        void reset() {
            first = true;
        }
    }

    public SubstitutionModel[] getSubstitutionModels() {
        return substitutionModels;
    }

    public final SubstitutionModel getSubstitutionModel() {
        if (substitutionModels.length != 1) {
            throw new IllegalStateException(
                    "Expected exactly one substitution model, but found "
                            + substitutionModels.length);
        }
        return substitutionModels[0];
    }

    public final int getSubstitutionModelCount() {
        return substitutionModels.length;
    }

    public final SubstitutionModel getSubstitutionModel(int index) {
        return substitutionModels[index];
    }

    public AbstractPopulationSizeModel getPopulationSizeModel() {
        return populationSizeModel;
    }

    @Override
    public final double getLogLikelihood() {
        onLikelihoodRequested();

        if (!likelihoodKnown) {
            long timingToken = beforeLikelihoodEvaluation();
            logLikelihood = calculateLogLikelihood();
            afterLikelihoodEvaluation(timingToken);
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    /**
     * Performs the engine-specific likelihood evaluation. Implementations
     * should return {@link Double#NEGATIVE_INFINITY} for an ordinary
     * numerically invalid proposal, rather than throwing or touching
     * {@link #likelihoodKnown} directly.
     */
    protected abstract double calculateLogLikelihood();

    /**
     * Lets a gradient or ancestral-state-reconstruction evaluation (which
     * computes the likelihood as a side effect of computing something else)
     * populate the shared likelihood cache, instead of touching {@link
     * #logLikelihood}/{@link #likelihoodKnown} directly.
     */
    protected final void cacheLogLikelihood(double value) {
        logLikelihood = value;
        likelihoodKnown = true;
    }

    protected final void invalidateLikelihood() {
        likelihoodKnown = false;
    }

    // Optional profiling hooks (see BastaLikelihood's COUNT_TOTAL_OPERATIONS
    // counters): a subclass that doesn't need them pays nothing extra, since
    // the default bodies are empty/zero and getLogLikelihood() is final so
    // every subclass's profiling goes through the same three call sites.
    protected void onLikelihoodRequested() {
    }

    protected long beforeLikelihoodEvaluation() {
        return 0L;
    }

    protected void afterLikelihoodEvaluation(long timingToken) {
    }

    @Override
    public final void makeDirty() {
        likelihoodKnown = false;
        invalidateDerivedState();
        makeDirtyInternal();
    }

    /** Invalidates engine-specific structural caches (event tape, transition matrices, etc.). */
    protected abstract void makeDirtyInternal();

    /**
     * Invalidates outputs derived from the likelihood but not part of its
     * own cache -- gradients, ancestral-state reconstructions/sensitivity
     * scores. Deliberately not a shared field: BASTA's reconstruction and
     * MASCOT's adjoint sensitivity scores are different quantities with
     * different validity semantics, so each subclass names and owns its own
     * cache flag(s); this hook is just the guaranteed call site for
     * invalidating them.
     */
    protected void invalidateDerivedState() {
    }

    @Override
    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        handleStructuredModelChangedEvent(model, object, index);
        likelihoodKnown = false;
        invalidateDerivedState();
        fireModelChanged();
    }

    protected abstract void handleStructuredModelChangedEvent(Model model, Object object, int index);

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        handleStructuredVariableChangedEvent(variable, index, type);
        likelihoodKnown = false;
        invalidateDerivedState();
        fireModelChanged();
    }

    protected void handleStructuredVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Unexpected variable change event: " + variable.getVariableName());
    }

    @Override
    protected final void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storeStructuredState();
    }

    /** Engine-specific store hook; the shared logLikelihood/likelihoodKnown snapshot above is not repeated here. */
    protected void storeStructuredState() {
    }

    @Override
    protected final void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        restoreStructuredState();
    }

    protected void restoreStructuredState() {
    }

    @Override
    protected void acceptState() {
        // Nothing to do: current fields already reflect the accepted state.
        // Both existing engines are no-ops here; a future one that isn't can override.
    }

    /** Matches both engines' previous default {@code getReport()} body; BASTA layers profiling/delegate output on top of it. */
    protected final String getDefaultReport() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    @Override
    public String getReport() {
        return getDefaultReport();
    }

    public double[] getGradientLogDensity(StructuredCoalescentLikelihoodGradient wrt) {
       throw new UnsupportedOperationException("Supported only in subclasses.");
    }
}
