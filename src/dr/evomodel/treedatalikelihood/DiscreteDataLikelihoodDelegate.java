package dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.*;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations.*;
import dr.evomodel.treedatalikelihood.preorder.DiscretePreOrderType;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.logging.Logger;

/**
 * @author Filippo Monti
 */
public class DiscreteDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate, Citable,
        PostOrderMessageProvider, PreOrderMessageProvider, GradientDataLikelihoodDelegate {

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel");
    private static final boolean COUNT_CALCULATIONS = true;
    private static final double DEFAULT_SCALING_FLOOR = 1.0e-200;
    private static final double DEFAULT_SCALING_CEILING = 1.0e200;

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        return TreeTraversal.TraversalType.POST_ORDER;
    }

    // -------------------------------------------------------------------------
    // Partial transforms used only for optional exported caches
    // -------------------------------------------------------------------------

    public interface PartialTransform {
        String getName();

        void transform(double[] input, double[] output);

        PartialTransform IDENTITY = new PartialTransform() {
            @Override
            public String getName() {
                return "identity";
            }

            @Override
            public void transform(double[] input, double[] output) {
                System.arraycopy(input, 0, output, 0, input.length);
            }
        };
    }

    // -------------------------------------------------------------------------
    // Cache settings
    // -------------------------------------------------------------------------

    public static final class CacheSettings {
        public final boolean usePreOrder;
        public final boolean cacheBranchStartPostOrder;
        public final boolean cacheBranchEndPostOrder;
        public final boolean cacheBranchStartPreOrder;
        public final boolean cacheBranchEndPreOrder;
        public final boolean cacheTransformedPostOrder;
        public final boolean cacheTransformedPreOrder;
        public final boolean applyPatternScaling;

        private CacheSettings(Builder b) {
            this.usePreOrder = b.usePreOrder;
            this.cacheBranchStartPostOrder = b.cacheBranchStartPostOrder;
            this.cacheBranchEndPostOrder = b.cacheBranchEndPostOrder;
            this.cacheBranchStartPreOrder = b.cacheBranchStartPreOrder;
            this.cacheBranchEndPreOrder = b.cacheBranchEndPreOrder;
            this.cacheTransformedPostOrder = b.cacheTransformedPostOrder;
            this.cacheTransformedPreOrder = b.cacheTransformedPreOrder;
            this.applyPatternScaling = b.applyPatternScaling;
        }

        public static Builder newBuilder() {
            return new Builder();
        }

        public static final class Builder {
            private boolean usePreOrder = false;
            private boolean cacheBranchStartPostOrder = false;
            private boolean cacheBranchEndPostOrder = true;
            private boolean cacheBranchStartPreOrder = false;
            private boolean cacheBranchEndPreOrder = false;
            private boolean cacheTransformedPostOrder = false;
            private boolean cacheTransformedPreOrder = false;
            private boolean applyPatternScaling = true;

            public Builder usePreOrder(boolean value) {
                this.usePreOrder = value;
                return this;
            }

            public Builder cacheBranchStartPostOrder(boolean value) {
                this.cacheBranchStartPostOrder = value;
                return this;
            }

            public Builder cacheBranchEndPostOrder(boolean value) {
                this.cacheBranchEndPostOrder = value;
                return this;
            }

            public Builder cacheBranchStartPreOrder(boolean value) {
                this.cacheBranchStartPreOrder = value;
                return this;
            }

            public Builder cacheBranchEndPreOrder(boolean value) {
                this.cacheBranchEndPreOrder = value;
                return this;
            }

            public Builder cacheTransformedPostOrder(boolean value) {
                this.cacheTransformedPostOrder = value;
                return this;
            }

            public Builder cacheTransformedPreOrder(boolean value) {
                this.cacheTransformedPreOrder = value;
                return this;
            }

            public Builder applyPatternScaling(boolean value) {
                this.applyPatternScaling = value;
                return this;
            }

            public CacheSettings build() {
                return new CacheSettings(this);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Core state
    // -------------------------------------------------------------------------

    private final Tree tree;
    private final PatternList patternList;
    private final DataType dataType;
    private final BranchModel branchModel;
    private final SiteRateModel siteRateModel;
    private final PostOrderRepresentation postOrderRepresentation;
    private final CacheSettings cacheSettings;
    private final PartialTransform postOrderTransform;
    private final PartialTransform preOrderTransform;
    private final boolean useAmbiguities;
    private final boolean preferGPU;
    private final PartialsRescalingScheme rescalingScheme;
    private final boolean delayRescalingUntilUnderflow;
    private final PreOrderSettings preOrderSettings;


    private final boolean preOrderEnabled;
    private boolean preOrderValid = false;
    private final PreOrderRepresentation preOrderRepresentation;
    private final DiscretePreOrderDelegate preOrderDelegate;
//    private final DiscreteLocalBranchUpdateEngine localBranchUpdateEngine;

    private final BidirectionalRepresentation bidirectionalRepresentation;

    private final Parameter siteAssignInd;
    private final Parameter partitionCat;

    private final int nodeCount;
    private final int tipCount;
    private final int patternCount;
    private final int categoryCount;
    private final int stateCount;

    private double[] patternWeights;
    private double[] storedPatternWeights;

    // [node][flattened(category,pattern,state)]
    private double[][] nodePartials;
    private double[][] storedNodePartials;

    // post-order scaling [node][pattern]
    private double[][] nodePatternLogScales;
    private double[][] storedNodePatternLogScales;

    // retained for API compatibility; not actively used in this first refactor
//    private double[][] nodePreOrderPatternLogScales;
//    private double[][] storedNodePreOrderPatternLogScales;

    // Optional caches exposed in STANDARD basis
    private final double[][] postOrderAtBranchStart;
    private final double[][] postOrderAtBranchEnd;
    private final double[][] preOrderAtBranchStart;
    private final double[][] preOrderAtBranchEnd;
    private final double[][] transformedPostOrder;
    private final double[][] transformedPreOrder;

    private final double[] branchLengths;
    private final int[] branchUpdateIndices;

    private final boolean[] nodePartialKnown;
    private final boolean[] nodePreOrderKnown;
    private final boolean[] postOrderStartKnown;
    private final boolean[] postOrderEndKnown;
    private final boolean[] preOrderStartKnown;
    private final boolean[] preOrderEndKnown;
    private final boolean[] transformedPostOrderKnown;
    private final boolean[] transformedPreOrderKnown;

    private final double[] patternLogLikelihoods;

    // Scratch vectors reused everywhere
    private final double[] leftPropagated;
    private final double[] rightPropagated;
    private final double[] tmpVectorA;
    private final double[] tmpVectorB;

    // Scratch buffer for exporting a whole node buffer to standard basis
    private final double[] tmpNodeStandardBuffer;

    private boolean updateSubstitutionModel = true;
    private boolean updateSiteModel = true;
    private boolean updateRootFrequency = true;
    private boolean computePostOrderStatisticsOnly = false;

    private long totalEvaluationCount = 0L;
    private long totalNodeOperationCount = 0L;
    private long totalBranchOperationCount = 0L;






    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------

    public DiscreteDataLikelihoodDelegate(Tree tree,
                                          PatternList patternList,
                                          BranchModel branchModel,
                                          SiteRateModel siteRateModel,
                                          boolean useAmbiguities,
                                          boolean preferGPU,
                                          PartialsRescalingScheme rescalingScheme,
                                          boolean delayRescalingUntilUnderflow,
                                          PreOrderSettings preOrderSettings,
                                          PostOrderRepresentation postOrderRepresentation,
                                          PartialTransform postOrderTransform, // TODO is this really useful?
                                          PartialTransform preOrderTransform, // TODO is this really useful?
                                          Parameter siteAssignInd,
                                          Parameter partitionCat) {
        super("DiscreteDataLikelihoodDelegate");

        this.tree = Objects.requireNonNull(tree, "tree");
        this.patternList = Objects.requireNonNull(patternList, "patternList");
        this.branchModel = Objects.requireNonNull(branchModel, "branchModel");
        this.siteRateModel = Objects.requireNonNull(siteRateModel, "siteRateModel");
//        this.postOrderRepresentation = Objects.requireNonNull(postOrderRepresentation, "postOrderRepresentation");
        this.postOrderTransform = postOrderTransform == null ? PartialTransform.IDENTITY : postOrderTransform;
        this.preOrderTransform = preOrderTransform == null ? PartialTransform.IDENTITY : preOrderTransform;
        this.useAmbiguities = useAmbiguities || patternList instanceof UncertainSiteList;
        this.preferGPU = preferGPU;
        this.rescalingScheme = rescalingScheme == null ? PartialsRescalingScheme.NONE : rescalingScheme;
        this.delayRescalingUntilUnderflow = delayRescalingUntilUnderflow;
        this.preOrderSettings = preOrderSettings == null
                ? new PreOrderSettings(false, false, false, false)
                : preOrderSettings;
        this.dataType = patternList.getDataType();
        this.patternCount = patternList.getPatternCount();
        this.stateCount = dataType.getStateCount();
        this.categoryCount = siteRateModel.getCategoryCount();
        this.nodeCount = tree.getNodeCount();
        this.tipCount = tree.getExternalNodeCount();
//        this.internalNodeCount = nodeCount - tipCount;
        this.siteAssignInd = siteAssignInd;
        this.partitionCat = partitionCat;

        this.branchLengths = new double[nodeCount];
        this.branchUpdateIndices = new int[nodeCount];


//        this.postOrderRepresentation = new StandardPostOrderRepresentation(this.postOrderBranchExecution, stateCount);
        this.postOrderRepresentation = Objects.requireNonNull(postOrderRepresentation, "postOrderRepresentation");

        if (this.postOrderRepresentation instanceof BidirectionalRepresentation) {
            this.bidirectionalRepresentation = (BidirectionalRepresentation) this.postOrderRepresentation;
        } else {
            this.bidirectionalRepresentation = null;
        }

        if (this.postOrderRepresentation == null) {
            throw new IllegalArgumentException("postOrderRepresentation must not be null");
        }
        if (this.postOrderRepresentation.getStateCount() != stateCount) {
            throw new IllegalArgumentException(
                    "Representation state count (" + this.postOrderRepresentation.getStateCount() +
                            ") does not match pattern state count (" + stateCount + ")");
        }

//        if (cacheSettings.usePreOrder ||
//                cacheSettings.cacheBranchStartPreOrder ||
//                cacheSettings.cacheBranchEndPreOrder ||
//                cacheSettings.cacheTransformedPreOrder) {
//            throw new UnsupportedOperationException(
//                    "Pre-order caching is intentionally not implemented in this first post-order generalization.");
//        }
        this.cacheSettings = buildDefaultCacheSettings(preOrderSettings.usePreOrder);

        addModel(this.branchModel);
        addModel(this.siteRateModel);
        if (siteAssignInd != null) {
            addVariable(siteAssignInd);
        }

        final int frequencyCount = branchModel.getRootFrequencyModel().getFrequencyCount();
        if (frequencyCount != stateCount) {
            throw new IllegalArgumentException(
                    "Pattern state count (" + stateCount + ") does not match root frequency count (" + frequencyCount + ")");
        }

        this.patternWeights = initialisePatternWeights();
        this.storedPatternWeights = new double[patternWeights.length];

        final int nodeBufferLength = flattenedPartialLength();
        this.nodePartials = new double[nodeCount][nodeBufferLength];
        this.storedNodePartials = new double[nodeCount][nodeBufferLength];

        this.nodePatternLogScales = new double[nodeCount][patternCount];
        this.storedNodePatternLogScales = new double[nodeCount][patternCount];

//        this.nodePreOrderPatternLogScales = new double[nodeCount][patternCount];
//        this.storedNodePreOrderPatternLogScales = new double[nodeCount][patternCount];

        this.postOrderAtBranchStart = cacheSettings.cacheBranchStartPostOrder ? new double[nodeCount][nodeBufferLength] : null;
        this.postOrderAtBranchEnd = cacheSettings.cacheBranchEndPostOrder ? new double[nodeCount][nodeBufferLength] : null;

        this.transformedPostOrder = cacheSettings.cacheTransformedPostOrder ? new double[nodeCount][nodeBufferLength] : null;

        this.preOrderAtBranchStart = cacheSettings.cacheBranchStartPreOrder ? new double[nodeCount][nodeBufferLength] : null;
        this.preOrderAtBranchEnd = cacheSettings.cacheBranchEndPreOrder ? new double[nodeCount][nodeBufferLength] : null;
        this.transformedPreOrder = cacheSettings.cacheTransformedPreOrder ? new double[nodeCount][nodeBufferLength] : null;

        this.preOrderStartKnown = preOrderAtBranchStart != null ? new boolean[nodeCount] : null;
        this.preOrderEndKnown = preOrderAtBranchEnd != null ? new boolean[nodeCount] : null;
        this.transformedPreOrderKnown = transformedPreOrder != null ? new boolean[nodeCount] : null;

        this.nodePartialKnown = new boolean[nodeCount];
        this.nodePreOrderKnown = new boolean[nodeCount];
        this.postOrderStartKnown = postOrderAtBranchStart != null ? new boolean[nodeCount] : null;
        this.postOrderEndKnown = postOrderAtBranchEnd != null ? new boolean[nodeCount] : null;

        this.transformedPostOrderKnown = transformedPostOrder != null ? new boolean[nodeCount] : null;

        this.patternLogLikelihoods = new double[patternCount];

        this.leftPropagated = new double[stateCount];
        this.rightPropagated = new double[stateCount];
        this.tmpVectorA = new double[stateCount];
        this.tmpVectorB = new double[stateCount];
        this.tmpNodeStandardBuffer = new double[nodeBufferLength];

        initialiseTipPartials();

        this.preOrderEnabled = cacheSettings.usePreOrder;
        if (preOrderEnabled) {
            if (bidirectionalRepresentation != null) {
                this.preOrderRepresentation = bidirectionalRepresentation;
            } else {
                throw new UnsupportedOperationException(
                        "Pre-order is enabled, but the supplied representation does not implement "
                                + "PreOrderRepresentation and no legacy pre-order branch execution was provided.");
            }

            this.preOrderDelegate =
                    new DiscretePreOrderDelegate(
                            tree,
                            branchLengths,
                            preOrderRepresentation,
                            this,
                            patternCount,
                            categoryCount,
                            cacheSettings.cacheBranchStartPreOrder,
                            cacheSettings.cacheBranchEndPreOrder
                    );

        } else {
            this.preOrderRepresentation = null;
            this.preOrderDelegate = null;
        }
//        this.localBranchUpdateEngine = null; //TODO
        LOGGER.info("\nCreating DiscreteDataLikelihoodDelegate");
        LOGGER.info("    Using Java discrete-state pruning engine.");
        LOGGER.info("    Post-order representation: " + this.postOrderRepresentation.getName());
        LOGGER.info("    Pre-order enabled: " + preOrderEnabled);
        LOGGER.info("    Post-order transform cache: " + cacheSettings.cacheTransformedPostOrder + " (" + this.postOrderTransform.getName() + ")");
        LOGGER.info("    Pattern scaling: " + cacheSettings.applyPatternScaling);
    }

    private static CacheSettings buildDefaultCacheSettings(boolean enablePreOrder) {
        return CacheSettings.newBuilder()
                .usePreOrder(enablePreOrder)
                .cacheBranchStartPostOrder(true)
                .cacheBranchEndPostOrder(true)
                .cacheBranchStartPreOrder(enablePreOrder)
                .cacheBranchEndPreOrder(enablePreOrder)
                .cacheTransformedPostOrder(false)
                .cacheTransformedPreOrder(false)
                .applyPatternScaling(true)
                .build();
    }

    public double getEffectiveBranchLength(int childNodeNumber) {
        return branchLengths[childNodeNumber];
    }

    public double[] getPatternWeights() {
        return patternWeights;
    }

    public double[] getCategoryWeights() {
        return siteRateModel.getCategoryProportions();
    }

    public BranchModel getBranchModel() {
        return branchModel;
    }

    public final SiteRateModel getSiteRateModel(){
        return this.siteRateModel;
    }

    // -------------------------------------------------------------------------
    // DataLikelihoodDelegate implementation
    // -------------------------------------------------------------------------

    @Override
    public void makeDirty() {
        updateSubstitutionModel = true;
        updateSiteModel = true;
        updateRootFrequency = true;
        invalidateAllCaches();
        postOrderRepresentation.markDirty();
        fireModelChanged();
        preOrderValid = false;
    }

    @Override
    public void storeState() {
        System.arraycopy(patternWeights, 0, storedPatternWeights, 0, patternWeights.length);
        copy2D(nodePartials, storedNodePartials);
        copy2D(nodePatternLogScales, storedNodePatternLogScales);
//        copy2D(nodePreOrderPatternLogScales, storedNodePreOrderPatternLogScales);
        postOrderRepresentation.storeState();
        if (preOrderDelegate != null) {
            preOrderDelegate.storeState();
        }
    }

    @Override
    public void restoreState() {
        double[] tmpWeights = patternWeights;
        patternWeights = storedPatternWeights;
        storedPatternWeights = tmpWeights;

        double[][] tmpNodePartials = nodePartials;
        nodePartials = storedNodePartials;
        storedNodePartials = tmpNodePartials;

        double[][] tmpNodePatternLogScales = nodePatternLogScales;
        nodePatternLogScales = storedNodePatternLogScales;
        storedNodePatternLogScales = tmpNodePatternLogScales;

//        double[][] tmpNodePreOrderPatternLogScales = nodePreOrderPatternLogScales;
//        nodePreOrderPatternLogScales = storedNodePreOrderPatternLogScales;
//        storedNodePreOrderPatternLogScales = tmpNodePreOrderPatternLogScales;

        updateSiteModel = true;
        updateRootFrequency = true;
        updateSubstitutionModel = true;
        invalidateAllCaches();
        initialiseTipPartials();
        postOrderRepresentation.restoreState();
        preOrderValid = false;
        if (preOrderDelegate != null) {
            preOrderDelegate.restoreState();
        }
    }
    public void invalidatePreOrderOnlyForDebug() {
        preOrderValid = false;
        Arrays.fill(nodePreOrderKnown, false);
        if (preOrderStartKnown != null) Arrays.fill(preOrderStartKnown, false);
        if (preOrderEndKnown != null) Arrays.fill(preOrderEndKnown, false);
        if (transformedPreOrderKnown != null) Arrays.fill(transformedPreOrderKnown, false);
//        for (int i = 0; i < nodeCount; i++) {
//            Arrays.fill(nodePreOrderPatternLogScales[i], 0.0);
//        }
        if (preOrderDelegate != null) {
            preOrderDelegate.makeDirty();
        }
    }

    @Override
    public double[] getSiteLogLikelihoods() {
        return Arrays.copyOf(patternLogLikelihoods, patternLogLikelihoods.length);
    }


    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations,
                                      List<NodeOperation> nodeOperations,
                                      int rootNodeNumber) throws LikelihoodException {

        if (siteAssignInd != null) {
            refreshPatternWeightsFromSiteAssignments();
        }
        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            final int nodeNumber = op.getBranchNumber();
//            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber(); //TODO check this
            branchLengths[nodeNumber] = op.getBranchLength();
            branchUpdateCount ++;
        }

        if (updateSubstitutionModel || updateSiteModel || updateRootFrequency) {
            postOrderRepresentation.updateForLikelihood();
        }

        final double[] categoryRates = siteRateModel.getCategoryRates();
        if (categoryRates == null) {
            return Double.NEGATIVE_INFINITY;
        }

        final double[] categoryWeights = siteRateModel.getCategoryProportions();
        final double[] rootFrequencies = branchModel.getRootFrequencyModel().getFrequencies();

        if (rootFrequencies.length != stateCount) {
            throw new IllegalStateException("Root frequencies must have length equal to state count");
        }

        markDirtyFromNodeOperations(nodeOperations);

        final NodeRef root = tree.getNode(rootNodeNumber);
        ensurePostOrder(root, categoryRates);

        if (computePostOrderStatisticsOnly) {
            updateSubstitutionModel = false;
            updateSiteModel = false;
            updateRootFrequency = false;
            if (COUNT_CALCULATIONS) {
                totalEvaluationCount++;
            }
            return Double.NaN;
        }

        final double logL = computeRootLogLikelihood(root, rootFrequencies, categoryWeights);

        updateSubstitutionModel = false;
        updateSiteModel = false;
        updateRootFrequency = false;

        if (COUNT_CALCULATIONS) {
            totalEvaluationCount++;
            totalNodeOperationCount += nodeOperations == null ? 0 : nodeOperations.size();
            totalBranchOperationCount += branchOperations == null ? 0 : branchOperations.size();
        }

        if (Double.isNaN(logL) || Double.isInfinite(logL)) {
            throw new LikelihoodUnderflowException();
        }

        return logL;
    }

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getTraitDim() {
        return patternCount;
    }

    @Override
    public int getPartitionCat() {
        return partitionCat == null ? 0 : (int) partitionCat.getParameterValue(0);
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        return RateRescalingScheme.NONE;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
        // no-op
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            if (k + 3 > operations.length) {
                throw new IllegalArgumentException("Insufficient operations array length");
            }
            operations[k++] = op.getNodeNumber();
            operations[k++] = op.getLeftChild();
            operations[k++] = op.getRightChild();
        }
        return k;
    }

    @Override
    public void setComputePostOrderStatisticsOnly(boolean computePostOrderStatisticsOnly) {
        this.computePostOrderStatisticsOnly = computePostOrderStatisticsOnly;
    }

    @Override
    public boolean providesPostOrderStatisticsOnly() {
        return true;
    }

    @Override
    public PreOrderSettings getPreOrderSettings() {
        return preOrderSettings;
    }

    @Override
    public boolean getPreferGPU() {
        return preferGPU;
    }

    @Override
    public boolean getUseAmbiguities() {
        return useAmbiguities;
    }

    @Override
    public PartialsRescalingScheme getRescalingScheme() {
        return rescalingScheme;
    }

    @Override
    public boolean getDelayRescalingUntilUnderflow() {
        return delayRescalingUntilUnderflow;
    }

    @Override
    public long getTotalCalculationCount() {
        return totalNodeOperationCount;
    }

    public void updatePostOrdersFromTreeDataLikelihood(TreeDataLikelihood treeDataLikelihood) {
        treeDataLikelihood.calculatePostOrderStatistics();
    }
    // -------------------------------------------------------------------------
    // Report / cite
    // -------------------------------------------------------------------------

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Java discrete-state tree likelihood delegate";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.AYRES_2019_BEAGLE);
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append(getClass().getSimpleName()).append('\n');
        sb.append("  evaluations = ").append(totalEvaluationCount).append('\n');
        sb.append("  node operations = ").append(totalNodeOperationCount).append('\n');
        sb.append("  branch operations = ").append(totalBranchOperationCount).append('\n');
        sb.append("  patterns = ").append(patternCount).append('\n');
        sb.append("  categories = ").append(categoryCount).append('\n');
        sb.append("  states = ").append(stateCount).append('\n');
        sb.append("  post-order representation = ").append(postOrderRepresentation.getName()).append('\n');
        sb.append("  pre-order enabled = ").append(preOrderEnabled).append('\n');
        if (preOrderEnabled && preOrderRepresentation != null) {
            sb.append("  pre-order representation = ").append(preOrderRepresentation.getName()).append('\n');
        }
        return sb.toString();
    }
    public boolean hasPreOrder() {
        return preOrderEnabled;
    }

//    public DiscreteLocalBranchUpdateEngine getLocalBranchUpdateEngine() {
//        if (!preOrderEnabled || localBranchUpdateEngine == null) {
//            throw new IllegalStateException("Local branch update engine is not enabled.");
//        }
//        ensurePreOrderComputed();
//        return localBranchUpdateEngine;
//    }

    // -------------------------------------------------------------------------
    // AbstractModel callbacks
    // -------------------------------------------------------------------------

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == siteRateModel) {
            updateSiteModel = true;
        } else if (model == branchModel) {
            updateSubstitutionModel = true;
            updateRootFrequency = true;
        }
        invalidateAllCaches();
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == siteAssignInd) {
            refreshPatternWeightsFromSiteAssignments();
            invalidateAllCaches();
        }
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    // -------------------------------------------------------------------------
    // External cache accessors
    // -------------------------------------------------------------------------

    public double[] getPostOrderAtBranchEnd(int nodeNumber) {
        requireCache(postOrderAtBranchEnd, "postOrderAtBranchEnd");
        return Arrays.copyOf(postOrderAtBranchEnd[nodeNumber], postOrderAtBranchEnd[nodeNumber].length);
    }

    public double[] getPostOrderAtBranchStart(int nodeNumber) {
        requireCache(postOrderAtBranchStart, "postOrderAtBranchStart");
        return Arrays.copyOf(postOrderAtBranchStart[nodeNumber], postOrderAtBranchStart[nodeNumber].length);
    }

    public double[] getPreOrderAtBranchStart(int nodeNumber) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchStart, "preOrderAtBranchStart");
        return Arrays.copyOf(preOrderAtBranchStart[nodeNumber], preOrderAtBranchStart[nodeNumber].length);
    }

    public double[] getPreOrderAtBranchEnd(int nodeNumber) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchEnd, "preOrderAtBranchEnd");
        return Arrays.copyOf(preOrderAtBranchEnd[nodeNumber], preOrderAtBranchEnd[nodeNumber].length);
    }

    public void getPostOrderAtBranchEndInto(int nodeNumber, double[] dest) {
        requireCache(postOrderAtBranchEnd, "postOrderAtBranchEnd");
        if (dest.length != stateCount) {
            throw new IllegalArgumentException("Destination length must equal stateCount");
        }
        System.arraycopy(postOrderAtBranchEnd[nodeNumber], 0, dest, 0, stateCount);
    }

    public void getPreOrderAtBranchStartInto(int nodeNumber, double[] dest) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchStart, "preOrderAtBranchStart");
        if (dest.length != stateCount) {
            throw new IllegalArgumentException("Destination length must equal stateCount");
        }
        System.arraycopy(preOrderAtBranchStart[nodeNumber], 0, dest, 0, stateCount);
    }

    public void getPreOrderAtBranchEndInto(int nodeNumber, double[] dest) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchEnd, "preOrderAtBranchEnd");
        if (dest.length != stateCount) {
            throw new IllegalArgumentException("Destination length must equal stateCount");
        }
        System.arraycopy(preOrderAtBranchEnd[nodeNumber], 0, dest, 0, stateCount);
    }


    @Override
    public void getPreorderPartials(int node, DiscretePreOrderType type, double[] out) {

        assert out.length >= categoryCount * patternCount * stateCount;

        ensurePreOrderComputed();

        if (type == DiscretePreOrderType.BOTTOM) {
            requireCache(preOrderAtBranchEnd, "preOrderAtBranchEnd");
            System.arraycopy(preOrderAtBranchEnd[node], 0, out, 0,
                    categoryCount * patternCount * stateCount);

        } else if (type == DiscretePreOrderType.TOP) {
            requireCache(preOrderAtBranchStart, "preOrderAtBranchStart");
            System.arraycopy(preOrderAtBranchStart[node], 0, out, 0,
                    categoryCount * patternCount * stateCount);
        } else {
            throw new IllegalArgumentException("Partial type is not yet implemented");
        }
    }

    @Override
    public void getPreOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchStart, "preOrderAtBranchStart");

        if (outStandardPartial.length != stateCount) {
            throw new IllegalArgumentException("Output length must equal stateCount");
        }

        final int off = offset(category, pattern, 0);
        System.arraycopy(preOrderAtBranchStart[childNodeNumber], off, outStandardPartial, 0, stateCount);
    }

    @Override
    public void getPreOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial) {
        ensurePreOrderComputed();
        requireCache(preOrderAtBranchEnd, "preOrderAtBranchEnd");

        if (outStandardPartial.length != stateCount) {
            throw new IllegalArgumentException("Output length must equal stateCount");
        }

        final int off = offset(category, pattern, 0);
        System.arraycopy(preOrderAtBranchEnd[childNodeNumber], off, outStandardPartial, 0, stateCount);
    }

    public double[] getTransformedPostOrder(int nodeNumber) {
        requireCache(transformedPostOrder, "transformedPostOrder");
        return Arrays.copyOf(transformedPostOrder[nodeNumber], transformedPostOrder[nodeNumber].length);
    }

    public double[] getTransformedPreOrder(int nodeNumber) {
        ensurePreOrderComputed();
        requireCache(transformedPreOrder, "transformedPreOrder");
        return Arrays.copyOf(transformedPreOrder[nodeNumber], transformedPreOrder[nodeNumber].length);
    }

//    public double[] getPreOrderBranchScales(int nodeNumber) {
//        ensurePreOrderComputed();
//        return Arrays.copyOf(nodePreOrderPatternLogScales[nodeNumber], nodePreOrderPatternLogScales[nodeNumber].length);
//    }

    public double[] getPostOrderBranchScales(int nodeNumber) {
        return Arrays.copyOf(nodePatternLogScales[nodeNumber], nodePatternLogScales[nodeNumber].length);
    }

    // -------------------------------------------------------------------------
    // Core likelihood logic
    // -------------------------------------------------------------------------

    private void ensurePostOrder(NodeRef node, double[] categoryRates) {
        final int nodeNumber = node.getNumber();
        if (nodePartialKnown[nodeNumber]) {
            return;
        }

        if (tree.isExternal(node)) {
            nodePartialKnown[nodeNumber] = true;
            materializeOptionalPostOrderCaches(nodeNumber);
            return;
        }

        final NodeRef leftChild = tree.getChild(node, 0);
        final NodeRef rightChild = tree.getChild(node, 1);
        ensurePostOrder(leftChild, categoryRates);
        ensurePostOrder(rightChild, categoryRates);

        final int leftNumber = leftChild.getNumber();
        final int rightNumber = rightChild.getNumber();
        final double leftLength = branchLengths[leftNumber];
//                tree.getBranchLength(leftChild);
        final double rightLength = branchLengths[rightNumber];
//                tree.getBranchLength(rightChild);

        final double[] nodeBuffer = nodePartials[nodeNumber];
        final double[] leftBuffer = nodePartials[leftNumber];
        final double[] rightBuffer = nodePartials[rightNumber];
        final double[] nodeScales = nodePatternLogScales[nodeNumber];
        final double[] leftScales = nodePatternLogScales[leftNumber];
        final double[] rightScales = nodePatternLogScales[rightNumber];

        for (int c = 0; c < categoryCount; c++) {
            final double leftEffectiveLength = leftLength * categoryRates[c];
            final double rightEffectiveLength = rightLength * categoryRates[c];

            for (int p = 0; p < patternCount; p++) {
                final int childOffset = offset(c, p, 0);
                final int parentOffset = childOffset;

                sliceInto(leftBuffer, childOffset, tmpVectorA);
                postOrderRepresentation.propagateToBranchTop(leftNumber, leftEffectiveLength, tmpVectorA, leftPropagated);

                sliceInto(rightBuffer, childOffset, tmpVectorB);
                postOrderRepresentation.propagateToBranchTop(rightNumber, rightEffectiveLength, tmpVectorB, rightPropagated);

                nodeScales[p] = leftScales[p] + rightScales[p];

                postOrderRepresentation.combineBranchTopPartials(leftPropagated, rightPropagated, tmpVectorA);
                System.arraycopy(tmpVectorA, 0, nodeBuffer, parentOffset, stateCount);

                if (cacheSettings.applyPatternScaling) {
                    nodeScales[p] += normalizePatternSlice(nodeBuffer, parentOffset);
                }
//                if (postOrderAtBranchStart != null) {
//                    System.arraycopy(leftPropagated, 0, postOrderAtBranchStart[leftNumber], childOffset, stateCount);
//                    postOrderStartKnown[leftNumber] = true;
//
//                    System.arraycopy(rightPropagated, 0, postOrderAtBranchStart[rightNumber], childOffset, stateCount);
//                    postOrderStartKnown[rightNumber] = true;
//                }
                if (postOrderAtBranchStart != null) {
                    postOrderRepresentation.toStandard(leftPropagated, tmpVectorA);
                    System.arraycopy(tmpVectorA, 0, postOrderAtBranchStart[leftNumber], childOffset, stateCount);
                    postOrderStartKnown[leftNumber] = true;

                    postOrderRepresentation.toStandard(rightPropagated, tmpVectorB);
                    System.arraycopy(tmpVectorB, 0, postOrderAtBranchStart[rightNumber], childOffset, stateCount);
                    postOrderStartKnown[rightNumber] = true;
                }

//                if (postOrderAtBranchStart != null) {
//                    sliceInto(nodeBuffer, parentOffset, tmpVectorA);
//                    postOrderRepresentation.toStandard(tmpVectorA, tmpVectorB);
//                    System.arraycopy(tmpVectorB, 0, postOrderAtBranchStart[nodeNumber], parentOffset, stateCount);
//                    postOrderStartKnown[nodeNumber] = true;
//                }
            }
        }

        nodePartialKnown[nodeNumber] = true;
        materializeOptionalPostOrderCaches(nodeNumber);
    }

    private double computeRootLogLikelihood(NodeRef root, double[] rootFrequencies, double[] categoryWeights) {
        final int rootNumber = root.getNumber();
        final double[] rootPartials = nodePartials[rootNumber];
        final double[] rootScales = nodePatternLogScales[rootNumber];

        double totalLogLikelihood = 0.0;
        for (int p = 0; p < patternCount; p++) {
            double patternLikelihood = 0.0;

            for (int c = 0; c < categoryCount; c++) {
                final int off = offset(c, p, 0);
                sliceInto(rootPartials, off, tmpVectorA);
                patternLikelihood += categoryWeights[c] *
                        postOrderRepresentation.rootContribution(rootFrequencies, tmpVectorA);
            }

            if (patternLikelihood <= 0.0 || Double.isNaN(patternLikelihood) || Double.isInfinite(patternLikelihood)) {
                return Double.NEGATIVE_INFINITY;
            }

            final double logSite = Math.log(patternLikelihood) + rootScales[p];
            patternLogLikelihoods[p] = logSite;
            totalLogLikelihood += patternWeights[p] * logSite;
        }

        return totalLogLikelihood;
    }

    private void materializeOptionalPostOrderCaches(int nodeNumber) {
        final double[] source = nodePartials[nodeNumber];

        if (postOrderAtBranchEnd == null && transformedPostOrder == null) {
            return;
        }

        convertNodeBufferToStandard(source, tmpNodeStandardBuffer);

        if (postOrderAtBranchEnd != null) {
            System.arraycopy(tmpNodeStandardBuffer, 0, postOrderAtBranchEnd[nodeNumber], 0, source.length);
            postOrderEndKnown[nodeNumber] = true;
        }

        if (transformedPostOrder != null) {
            transformNodeBuffer(tmpNodeStandardBuffer, transformedPostOrder[nodeNumber], postOrderTransform);
            transformedPostOrderKnown[nodeNumber] = true;
        }
    }

    @Override
    public void getPostOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial) {
        requireCache(postOrderAtBranchStart, "postOrderAtBranchStart");
        final int off = offset(category, pattern, 0);
        System.arraycopy(postOrderAtBranchStart[childNodeNumber], off, outStandardPartial, 0, stateCount);
    }

    @Override
    public void getPostOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial) {
        requireCache(postOrderAtBranchEnd, "postOrderAtBranchEnd");
        final int off = offset(category, pattern, 0);
        System.arraycopy(postOrderAtBranchEnd[childNodeNumber], off, outStandardPartial, 0, stateCount);
    }

    @Override
    public void getPostOrderBranchScalesInto(int nodeNumber, double[] dest) {
        if (dest.length != patternCount) {
            throw new IllegalArgumentException("Destination length must equal patternCount");
        }
        System.arraycopy(nodePatternLogScales[nodeNumber], 0, dest, 0, patternCount);
    }
//    public getPreOrderBranchScalesInto(int nodeNumber, double[] dest) {
//        if (dest.length != patternCount) {
//            throw new IllegalArgumentException("Destination length must equal patternCount");
//        }
//        System.arraycopy(nodePreOrderPatternLogScales[nodeNumber], 0, dest, 0, patternCount);
//    }
    public void getPreOrderBranchScalesInto(int nodeNumber, double[] dest) {
        ensurePreOrderComputed();
        preOrderDelegate.getPreOrderBranchScalesInto(nodeNumber, dest);
    }

    @Override
    public int getStateCount() {
        return stateCount;
    }

    @Override
    public int getPatternCount() {
        return patternCount;
    }

    @Override
    public int getCategoryCount() {
        return categoryCount;
    }

    // -------------------------------------------------------------------------
    // Tip initialisation
    // -------------------------------------------------------------------------

    private void initialiseTipPartials() {
        for (int i = 0; i < tipCount; i++) {
            final String taxonId = tree.getTaxonId(i);
            final int sequenceIndex = patternList.getTaxonIndex(taxonId);
            if (sequenceIndex < 0) {
                throw new RuntimeException(
                        new TaxonList.MissingTaxonException("Taxon, " + taxonId + ", in tree, " + tree.getId()
                                + ", is not found in patternList, " + patternList.getId()));
            }
            initialiseTipNode(i, sequenceIndex);
        }
    }

    private void initialiseTipNode(int nodeNumber, int sequenceIndex) {
        final double[] partials = nodePartials[nodeNumber];
        Arrays.fill(partials, 0.0);
        Arrays.fill(nodePatternLogScales[nodeNumber], 0.0);
//        Arrays.fill(nodePreOrderPatternLogScales[nodeNumber], 0.0);

        final double[] tipStateBuffer = tmpVectorA;
        final double[] partialBuffer = tmpVectorB;

        for (int p = 0; p < patternCount; p++) {
            fillTipStateVector(sequenceIndex, p, tipStateBuffer);
            postOrderRepresentation.initializeTipPartial(tipStateBuffer, partialBuffer);

            for (int c = 0; c < categoryCount; c++) {
                final int off = offset(c, p, 0);
                System.arraycopy(partialBuffer, 0, partials, off, stateCount);
            }
        }

        nodePartialKnown[nodeNumber] = true;
        materializeOptionalPostOrderCaches(nodeNumber);
    }

    private void fillTipStateVector(int sequenceIndex, int patternIndex, double[] dest) {
        Arrays.fill(dest, 0.0);

        if (patternList instanceof UncertainSiteList) {
            ((UncertainSiteList) patternList).fillPartials(sequenceIndex, patternIndex, dest, 0);
            return;
        }

        if (patternList.areUncertain()) {
            final double[] prob = patternList.getUncertainPatternState(sequenceIndex, patternIndex);
            System.arraycopy(prob, 0, dest, 0, stateCount);
            return;
        }

        final int state = patternList.getPatternState(sequenceIndex, patternIndex);
        final boolean[] stateSet = dataType.getStateSet(state);
        for (int s = 0; s < stateCount; s++) {
            dest[s] = stateSet[s] ? 1.0 : 0.0;
        }
    }

    // -------------------------------------------------------------------------
    // Pre-order Machinery
    // -------------------------------------------------------------------------

    public void ensurePreOrderComputed() {
        if (!preOrderEnabled) {
            throw new IllegalStateException("Pre-order is not enabled for this delegate.");
        }

        if (preOrderValid) {
            return;
        }

        // Ensure likelihood/postorder is current first
        // safest choice: force current likelihood evaluation if needed
        // assuming TreeDataLikelihood has already driven calculateLikelihood before this is called,
        // you may only need the internal postorder to be up-to-date.
        final double[] categoryRates = siteRateModel.getCategoryRates();
        final double[] rootFrequencies = branchModel.getRootFrequencyModel().getFrequencies();
        final int rootNodeNumber = tree.getRoot().getNumber();

        if (updateSubstitutionModel || updateSiteModel || updateRootFrequency) { //TODO check this better
            postOrderRepresentation.updateForLikelihood();
            if (!(postOrderRepresentation instanceof BidirectionalRepresentation) && preOrderRepresentation != null) {
                preOrderRepresentation.updateForLikelihood();
            }
        }

        ensurePostOrder(tree.getRoot(), categoryRates);

//        preOrderDelegate.makeDirty();
        preOrderDelegate.ensurePreOrder(rootNodeNumber, categoryRates, rootFrequencies);
        materializeOptionalPreOrderCaches();

        preOrderValid = true;
    }

    private void materializeOptionalPreOrderCaches() {
        if (!preOrderEnabled) {
            return;
        }

        for (int nodeNumber = 0; nodeNumber < nodeCount; nodeNumber++) {
            if (preOrderAtBranchStart != null) {
//                final double[] src = preOrderDelegate.getPreOrderAtBranchStart(nodeNumber);
//                System.arraycopy(src, 0, preOrderAtBranchStart[nodeNumber], 0, src.length);
                preOrderDelegate.getPreOrderAtBranchStartInto(nodeNumber, 0, 0, preOrderAtBranchStart[nodeNumber]);
                preOrderStartKnown[nodeNumber] = true;
            }

            if (preOrderAtBranchEnd != null) {
//                final double[] src = preOrderDelegate.getPreOrderAtBranchEnd(nodeNumber);
//                System.arraycopy(src, 0, preOrderAtBranchEnd[nodeNumber], 0, src.length);
                preOrderDelegate.getPreOrderAtBranchEndInto(nodeNumber, 0,0, preOrderAtBranchEnd[nodeNumber]);
                preOrderEndKnown[nodeNumber] = true;
            }

            if (transformedPreOrder != null) {
                final double[] sourceStandard =
                        preOrderAtBranchEnd != null
                                ? preOrderAtBranchEnd[nodeNumber]
                                : preOrderDelegate.getPreOrderAtBranchEnd(nodeNumber);

                transformNodeBuffer(sourceStandard, transformedPreOrder[nodeNumber], preOrderTransform);
                transformedPreOrderKnown[nodeNumber] = true;
            }
        }

//        for (int nodeNumber = 0; nodeNumber < nodeCount; nodeNumber++) {
//            nodePreOrderKnown[nodeNumber] = true;
//
//            System.arraycopy(
//                    preOrderDelegate.getPreOrderBranchScales(nodeNumber),
//                    0,
//                    nodePreOrderPatternLogScales[nodeNumber],
//                    0,
//                    patternCount
//            );
//        }
    }
    // -------------------------------------------------------------------------
    // Dirtying / pattern weights
    // -------------------------------------------------------------------------

    private void invalidateAllCaches() {
        preOrderValid = false;
        Arrays.fill(nodePartialKnown, false);
        Arrays.fill(nodePreOrderKnown, false);
        if (preOrderDelegate != null) {
            preOrderDelegate.makeDirty();
        }
        if (postOrderStartKnown != null) Arrays.fill(postOrderStartKnown, false);
        if (postOrderEndKnown != null) Arrays.fill(postOrderEndKnown, false);
        if (preOrderStartKnown != null) Arrays.fill(preOrderStartKnown, false);
        if (preOrderEndKnown != null) Arrays.fill(preOrderEndKnown, false);
        if (transformedPostOrderKnown != null) Arrays.fill(transformedPostOrderKnown, false);
        if (transformedPreOrderKnown != null) Arrays.fill(transformedPreOrderKnown, false);

//        for (int i = 0; i < nodeCount; i++) {
//            Arrays.fill(nodePreOrderPatternLogScales[i], 0.0);
//        }

        for (int i = 0; i < tipCount; i++) {
            nodePartialKnown[i] = true;
            materializeOptionalPostOrderCaches(i);
        }
    }

    private void markDirtyFromNodeOperations(List<NodeOperation> nodeOperations) {
        if (nodeOperations == null || nodeOperations.isEmpty()) {
            if (updateSubstitutionModel || updateSiteModel || updateRootFrequency) {
                invalidateAllCaches();
            }
            return;
        }

        for (NodeOperation op : nodeOperations) {
            invalidateNodeAndAncestors(op.getNodeNumber());
        }
    }

    private void invalidateNodeAndAncestors(int nodeNumber) {
        preOrderValid = false;
        NodeRef node = tree.getNode(nodeNumber);
        while (node != null) {
            final int n = node.getNumber();
            nodePartialKnown[n] = false;
            nodePreOrderKnown[n] = false;

            if (postOrderStartKnown != null) postOrderStartKnown[n] = false;
            if (postOrderEndKnown != null) postOrderEndKnown[n] = false;
            if (preOrderStartKnown != null) preOrderStartKnown[n] = false;
            if (preOrderEndKnown != null) preOrderEndKnown[n] = false;
            if (transformedPostOrderKnown != null) transformedPostOrderKnown[n] = false;
            if (transformedPreOrderKnown != null) transformedPreOrderKnown[n] = false;

            if (tree.isRoot(node)) {
                break;
            }
            node = tree.getParent(node);
        }

        for (int i = 0; i < tipCount; i++) {
            nodePartialKnown[i] = true;
        }
    }

    private double[] initialisePatternWeights() {
        if (siteAssignInd == null) {
            return Arrays.copyOf(patternList.getPatternWeights(), patternCount);
        }
        final double[] weights = new double[patternCount];
        refreshPatternWeightsFromSiteAssignments(weights);
        return weights;
    }

    private void refreshPatternWeightsFromSiteAssignments() {
        refreshPatternWeightsFromSiteAssignments(patternWeights);
    }

    private void refreshPatternWeightsFromSiteAssignments(double[] weights) {
        Arrays.fill(weights, 0.0);
        final int selectedPartition = getPartitionCat();
        for (int i = 0; i < siteAssignInd.getSize(); i++) {
            final int cat = (int) siteAssignInd.getParameterValue(i);
            if (cat == selectedPartition) {
                final int patternIndex = patternList.getPatternIndex(i);
                weights[patternIndex] += 1.0;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Buffer utilities
    // -------------------------------------------------------------------------

    private int flattenedPartialLength() {
        return categoryCount * patternCount * stateCount;
    }

    private int offset(int category, int pattern, int state) {
        return ((category * patternCount) + pattern) * stateCount + state;
    }

    private void sliceInto(double[] array, int off, double[] dest) {
        System.arraycopy(array, off, dest, 0, stateCount);
    }

    private double normalizePatternSlice(double[] buffer, int off) {
        double max = 0.0;
        for (int s = 0; s < stateCount; s++) {
            max = Math.max(max, Math.abs(buffer[off + s]));
        }

        if (max == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }

        if (max < DEFAULT_SCALING_FLOOR || max > DEFAULT_SCALING_CEILING) {
            for (int s = 0; s < stateCount; s++) {
                buffer[off + s] /= max;
            }
            return Math.log(max);
        }

        return 0.0;
    }

    private void transformNodeBuffer(double[] source, double[] target, PartialTransform transform) {
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                sliceInto(source, off, tmpVectorA);
                transform.transform(tmpVectorA, tmpVectorB);
                System.arraycopy(tmpVectorB, 0, target, off, stateCount);
            }
        }
    }

    private void convertNodeBufferToStandard(double[] source, double[] target) {
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                sliceInto(source, off, tmpVectorA);
                postOrderRepresentation.toStandard(tmpVectorA, tmpVectorB);
                System.arraycopy(tmpVectorB, 0, target, off, stateCount);
            }
        }
    }

    private static void copy2D(double[][] src, double[][] dst) {
        for (int i = 0; i < src.length; i++) {
            System.arraycopy(src[i], 0, dst[i], 0, src[i].length);
        }
    }

    private static void requireCache(Object cache, String name) {
        if (cache == null) {
            throw new IllegalStateException("Cache not enabled: " + name);
        }
    }
}
