/*
 * BeagleDataLikelihoodDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood;

import beagle.*;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchmodel.TransitionMatrixProviderBranchModel;
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel.EpochTimeProvider;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import java.util.Collections;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodel.treedatalikelihood.BeagleFunctionality.*;
import static dr.evomodel.treedatalikelihood.SubstitutionModelDelegate.BUFFER_POOL_SIZE_DEFAULT;

/**
 * BeagleDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate that uses BEAGLE
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 */

public class BeagleDataLikelihoodDelegate extends AbstractModel implements
        DataLikelihoodDelegate, DiscreteProcessAlongTree, TipStateAccessor, Citable, GradientDataLikelihoodDelegate {

    private static final boolean COUNT_CALCULATIONS = true; // keep a cumulative total of number of computations

    private static final boolean RESCALING_OFF = false; // a debugging switch

    private static final boolean DEBUG = false; // write debug information to stdOut

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    private static final String RESOURCE_AUTO_PROPERTY = "beagle.resource.auto";
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    private static final String DELAY_SCALING_PROPERTY = "beagle.delay.scaling";
    private static final String EXTRA_BUFFER_COUNT_PROPERTY = "beagle.extra.buffer.count";
    private static final String FORCE_VECTORIZATION = "beagle.force.vectorization";
    private static final String THREAD_COUNT = "beagle.thread.count";
    // keep BEAGLE buffers for every augmented epoch node the tree can have, instead of growing them on demand
    private static final String PREALLOCATE_AUGMENTED_PROPERTY = "beagle.augmented.preallocate";
    // set to false to compute nodes with one ghost child through an all-unknown tip instead of a one-child operation
    public static final String GHOST_ONE_CHILD_PROPERTY = "beagle.ghost.one.child";
    // set to false to compute nodes with one ghost child instead of joining the branches through them
    public static final String GHOST_MERGE_PROPERTY = "beagle.ghost.merge";
    // with useSpectralRepresentation, ask BEAGLE to reuse post-order messages and pre-order push-downs
    private static final String SPECTRAL_CACHE_PROPERTY = "beagle.spectral.cache";
    private static final long SPECTRAL_CACHE_FLAG = 1L << 43; // BEAGLE_FLAG_SPECTRAL_CACHE; BeagleFlag.SPECTRAL_CACHE
                                                              // once lib/beagle.jar has it
    // on CPU, ask BEAGLE to compute the adjoint derivative for all branches at once
    private static final String ADJOINT_BATCHED_PROPERTY = "beagle.adjoint.batched";
    private static final long ADJOINT_BATCHED_FLAG = 1L << 42; // BEAGLE_FLAG_ADJOINT_BATCHED; BeagleFlag.ADJOINT_BATCHED
                                                                // once lib/beagle.jar has it

    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

    public static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Long> preferredOrder = null;
    private static List<Long> requiredOrder = null;
    private static List<String> scalingOrder = null;
    private static List<Integer> extraBufferOrder = null;

    // Default frequency for complete recomputation of scaling factors under the 'dynamic' scheme
    private static final int RESCALE_FREQUENCY = 10000;
    private static final int RESCALE_TIMES = 1;

    // count the number of partial likelihood and matrix updates
    private long totalMatrixUpdateCount = 0;
    private long totalPartialsUpdateCount = 0;
    private long totalEvaluationCount = 0;
    private boolean releaseSingleton = true;

    /**
     *
     * @param tree Used for configuration - shouldn't be watched for changes
     * @param branchModel Specifies substitution model for each branch
     * @param patternList List of patterns
     * @param siteRateModel Specifies rates per site
     * @param useAmbiguities Whether to respect state ambiguities in data
     * @param settings pre-order settings
     */
    public BeagleDataLikelihoodDelegate(Tree tree,
                                        PatternList patternList,
                                        BranchModel branchModel,
                                        SiteRateModel siteRateModel,
                                        boolean useAmbiguities,
                                        boolean preferGPU,
                                        PartialsRescalingScheme rescalingScheme,
                                        boolean delayRescalingUntilUnderflow,
                                        PreOrderSettings settings,
                                        Parameter siteAssignInd,
                                        Parameter partitionCat) {
        this(tree, patternList, branchModel, siteRateModel, useAmbiguities, preferGPU, rescalingScheme,
                delayRescalingUntilUnderflow, settings, siteAssignInd, partitionCat, false);
    }

    /**
     * @param allowGhostTaxa if true, taxa of the tree that are not in the pattern list are ghosts: they have no data
     *                       in this partition, so nothing is computed for any subtree made only of them. If false,
     *                       such a taxon is an error.
     */
    public BeagleDataLikelihoodDelegate(Tree tree,
                                        PatternList patternList,
                                        BranchModel branchModel,
                                        SiteRateModel siteRateModel,
                                        boolean useAmbiguities,
                                        boolean preferGPU,
                                        PartialsRescalingScheme rescalingScheme,
                                        boolean delayRescalingUntilUnderflow,
                                        PreOrderSettings settings,
                                        Parameter siteAssignInd,
                                        Parameter partitionCat,
                                        boolean allowGhostTaxa) {

        super("BeagleDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nCreating BEAGLE DataLikelihood Delegate");
        setId(patternList.getId());

        this.dataType = patternList.getDataType();
        this.patternList = patternList;
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

        if(siteAssignInd != null) {
            this.partitionCat = partitionCat;
            this.siteAssignInd = siteAssignInd;
            addVariable(siteAssignInd);
        }

        // Check for matching state counts
        int stateCount2 = branchModel.getRootFrequencyModel().getFrequencyCount();
        if (stateCount != stateCount2) {
            throw new IllegalArgumentException("Pattern state count (" + stateCount
                    + ") does not match substitution model state count (" + stateCount2 + ")");
        }

        if(siteAssignInd == null) {
            patternWeights = patternList.getPatternWeights();
        }else{
            patternWeights = new double[patternCount];
            // following line commented out because patternWeights array is initialized with default values of 0.0
            //Arrays.fill(patternWeights, 0.0);
            for(int i = 0; i < siteAssignInd.getSize(); i++){
                int cat = (int) siteAssignInd.getParameterValue(i);
                if((int) partitionCat.getParameterValue(0) == cat) {
                    int pwIndex = patternList.getPatternIndex(i);
                    patternWeights[pwIndex] = patternWeights[pwIndex] + 1;
                }
            }
            storedPatternWeights = new double[patternWeights.length];
        }

        this.branchModel = branchModel;
        if (!settings.useAugmentedEpochNodes) { // see below
            addModel(this.branchModel);
        }

        this.siteRateModel = siteRateModel;
        addModel(this.siteRateModel);

        this.categoryCount = this.siteRateModel.getCategoryCount();

        nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        isGhostTip = allowGhostTaxa ? findGhostTips(tree, patternList) : null;
        if (isGhostTip != null) {
            final int ghostCount = countGhosts(isGhostTip);
            if (ghostCount == tipCount) {
                throw new IllegalArgumentException("The pattern list, " + patternList.getId() +
                        ", has none of the taxa in the tree, " + tree.getId());
            }
            final String problem = ghostProblem(settings);
            if (problem != null) {
                throw new IllegalArgumentException("The tree, " + tree.getId() + ", has " + ghostCount +
                        " taxa that are not in the pattern list, " + patternList.getId() + ", but " + problem);
            }
        }

        // BEAGLE's spectral implementations keep the eigen information of each branch instead of its transition
        // matrix, so there are no matrices to convolve
        if (settings.useSpectralRepresentation && branchModel.requiresMatrixConvolution() &&
                !settings.useAugmentedEpochNodes) {
            throw new IllegalArgumentException("The spectral representation cannot convolve transition matrices; " +
                    "epoch models need useAugmentedEpochNodes=\"true\" with useSpectralRepresentation=\"true\"");
        }
        if (settings.useSpectralRepresentation) {
            checkMatrixVectorSupported("The spectral representation (useSpectralRepresentation=\"true\")");
        }

        // With augmented epoch nodes, degree-2 nodes lie along the branches and need their own buffers
        if (settings.useAugmentedEpochNodes) {
            if (!(branchModel instanceof EpochBranchModel)) {
                throw new IllegalArgumentException("Augmented epoch nodes require an epoch branch model");
            }
            if (settings.branchInfinitesimalDerivative) {
                throw new IllegalArgumentException("Augmented epoch nodes do not yet support " +
                        "branch infinitesimal derivatives");
            }
            checkMatrixVectorSupported("Augmented epoch nodes (useAugmentedEpochNodes=\"true\")");

            // An epoch branch model passes on every change of the tree, because the branches that a substitution
            // model applies to depend on the node heights, and then everything is updated. With degree-2 nodes the
            // heights only change the branch lengths, which the likelihood already updates where they change.
            // So, watch the parts of the branch model that matter instead of the branch model.
            for (SubstitutionModel substitutionModel : branchModel.getSubstitutionModels()) {
                addModel(substitutionModel);
            }
            addModel(branchModel.getRootFrequencyModel());

            // registered as a model so that it is restored when a proposal is rejected
            epochTimeProvider = new EpochTimeProvider.ParameterWrapper(
                    ((EpochBranchModel) branchModel).getEpochTimes());
            addModel(epochTimeProvider);

            epochDelegate = new EpochSubstitutionModelDelegate(tree, branchModel, epochTimeProvider, settings);
        } else {
            epochDelegate = null;
        }
        registry = (epochDelegate == null) ? null : epochDelegate.getAugmentedNodeRegistry();
        growable = epochDelegate != null && !Boolean.getBoolean(PREALLOCATE_AUGMENTED_PROPERTY);

        final int augmentedNodeCapacity = (epochDelegate == null) ? 0 :
                epochDelegate.getAugmentedNodeRegistry().getCapacity();

        branchUpdateIndices = new int[nodeCount + augmentedNodeCapacity];
        branchLengths = new double[nodeCount + augmentedNodeCapacity];
        scaleBufferIndices = new int[internalNodeCount];
        storedScaleBufferIndices = new int[internalNodeCount];

        operations = new int[(internalNodeCount + augmentedNodeCapacity) * Beagle.OPERATION_TUPLE_SIZE];

        firstRescaleAttempt = true;

        isRestored = false;

        this.settings = settings;

        try {

            // if we are using ambiguities then we don't use tip partials
            compactPartialsCount = useAmbiguities ? 0 : tipCount;

            // one partials buffer for each tip and two for each internal node (for store restore)
            if (epochDelegate != null) {
                // then, last so that allocating more augmented nodes appends, a block for each augmented node: its
                // two partials and, with pre-order, its pre-order partial (those of the original nodes come first)
                augmentedPartialHelper = new AugmentedBufferIndexHelper(nodeCount, tipCount, augmentedNodeCapacity,
                        2 * nodeCount - tipCount + (settings.usePreOrder ? nodeCount : 0),
                        settings.usePreOrder ? 3 : 2);
                partialBufferHelper = augmentedPartialHelper;
            } else {
                augmentedPartialHelper = null;
                partialBufferHelper = new BufferIndexHelper(nodeCount + augmentedNodeCapacity, tipCount);
            }

            // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
            scaleBufferHelper = new BufferIndexHelper(getSingleScaleBufferCount(), 0);

            if (extraBufferOrder == null) {
                extraBufferOrder = parseSystemPropertyIntegerArray(EXTRA_BUFFER_COUNT_PROPERTY);
            }

            int extraBufferCount = BUFFER_POOL_SIZE_DEFAULT;
            if (extraBufferOrder.size() > 0) {
                extraBufferCount = extraBufferOrder.get(instanceCount % extraBufferOrder.size());
            }

            if (epochDelegate != null) {
                evolutionaryProcessDelegate = epochDelegate;
            } else if (settings.useRewardAwareBranchModelDelegate) {
                evolutionaryProcessDelegate = new RewardAwareSubstitutionModelDelegate(tree,
                        (TransitionMatrixProviderBranchModel) branchModel, 0, extraBufferCount, settings);
            }  else if (settings.branchInfinitesimalDerivative) {
                    evolutionaryProcessDelegate = new SubstitutionModelDelegate(tree, branchModel, 0,
                            extraBufferCount, settings);
            } else {
                    if (branchModel.getSubstitutionModels().size() == 1) {
                        evolutionaryProcessDelegate = new HomogenousSubstitutionModelDelegate(tree, branchModel);
                    } else {
                        // use a more general delegate that allows different substitution models on different branches and
                        // can do matrix convolution.

                        // TODO: the constructor should take the delegate and the delegate should wrap the branchModel
                        evolutionaryProcessDelegate = new SubstitutionModelDelegate(tree, branchModel, 0,
                                extraBufferCount, settings);
                    }
            }

            int numPartials = partialBufferHelper.getBufferCount();
            int numScaleBuffers = scaleBufferHelper.getBufferCount();
            int numMatrices = evolutionaryProcessDelegate.getMatrixBufferCount();

            if (epochDelegate != null) {
                // the augmented layouts already include the pre-order partials, the gradient scale buffers and the
                // cached matrices; when they can grow, there are no buffers for augmented nodes until they are used
                final int allocated = growable ? 0 : augmentedNodeCapacity;
                numPartials = partialsFor(allocated);
                numMatrices = matricesFor(allocated);
                numScaleBuffers = scaleFor(allocated);
                setAugmentedAllocation(allocated);

            // one partial buffer for root node and two for each node including tip nodes (for store restore)
            } else if (settings.usePreOrder){
                // (pre-order buffers are also needed for the augmented nodes)
                numPartials += nodeCount + augmentedNodeCapacity;
                // Reserve one more scale buffer per node (indexed directly by node
                // number, not double-buffered/store-restore -- rebuilt fresh every
                // gradient evaluation), holding the negated cumulative sum of every
                // ANCESTOR's own individual rescaling. Used by
                // AbstractBeagleGradientDelegate/SpectralBeagleCrossProductDelegate
                // to correct the adjoint cross-product gradient for rescaling -- see
                // beagle.h's beagleCalculateAdjointCrossProductDerivative and
                // AbstractBeagleGradientDelegate.buildNegAncestorCumulativeScaleBuffer
                // for why this single term suffices (no separate post-order-cumulative
                // or root-cumulative buffer needed). Root's own slot goes unused; the
                // waste is one scale buffer's worth of memory (O(patternCount)
                // doubles) and not worth avoiding here.
                numScaleBuffers += nodeCount + augmentedNodeCapacity;
                numMatrices += evolutionaryProcessDelegate.getCachedMatrixBufferCount(settings);
            }

            // Attempt to get the resource order from the System Property
            if (resourceOrder == null) {
                resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
            }
            if (preferredOrder == null) {
                preferredOrder = parseSystemPropertyLongArray(PREFERRED_FLAGS_PROPERTY);
            }
            if (requiredOrder == null) {
                requiredOrder = parseSystemPropertyLongArray(REQUIRED_FLAGS_PROPERTY);
            }
            if (scalingOrder == null) {
                scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
            }

            // first set the rescaling scheme to use from the parser
            this.rescalingScheme = rescalingScheme;
            this.delayRescalingUntilUnderflow = delayRescalingUntilUnderflow;
            this.prefersGPU = preferGPU;
            this.useAmbiguities = useAmbiguities;

            int[] resourceList = null;
            long preferenceFlags = 0;
            long requirementFlags = 0;

            if (scalingOrder.size() > 0) {
                this.rescalingScheme = PartialsRescalingScheme.parseFromString(
                        scalingOrder.get(instanceCount % scalingOrder.size()));
            }

            if (resourceOrder.size() > 0) {
                // added the zero on the end so that a CPU is selected if requested resource fails
                resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
                if (!IS_CPU_RESOURCE(resourceList[0])) {
                    preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
                }
            }
            if (preferGPU) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }

            if (preferredOrder.size() > 0) {
                preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
            }

            if (requiredOrder.size() > 0) {
                requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
            }

            if (settings.useSpectralRepresentation) {
                requirementFlags |= BeagleFlag.SPECTRAL_REPRESENTATION.getMask();
            }

            if (epochDelegate != null) {
                // degree-2 nodes have Beagle.NONE as the second child, which only the CPU implementations accept
                requirementFlags |= BeagleFlag.FRAMEWORK_CPU.getMask();
            }


            // Define default behaviour here
            if (this.rescalingScheme == PartialsRescalingScheme.DEFAULT) {
                //if GPU: the default is dynamic scaling in BEAST
                if (resourceList != null && resourceList[0] > 1) {
                    this.rescalingScheme = DEFAULT_RESCALING_SCHEME;
                } else { // if CPU: just run as fast as possible
//                    this.rescalingScheme = PartialsRescalingScheme.NONE;
                    // Dynamic should run as fast as none until first underflow
                    this.rescalingScheme = DEFAULT_RESCALING_SCHEME;
                }
            }

            // to keep behaviour of the delayed scheme (always + delay)...
            if (this.rescalingScheme == PartialsRescalingScheme.DELAYED) {
                this.delayRescalingUntilUnderflow = true;
                this.rescalingScheme = PartialsRescalingScheme.ALWAYS;
            }

            if (this.rescalingScheme == PartialsRescalingScheme.AUTO) {
                if (epochDelegate != null) {
                    throw new IllegalArgumentException("Augmented epoch nodes do not support auto rescaling");
                }
                preferenceFlags |= BeagleFlag.SCALING_AUTO.getMask();
                useAutoScaling = true;
            } else {
//                preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
            }

            String r = System.getProperty(RESCALE_FREQUENCY_PROPERTY);
            if (r != null) {
                rescalingFrequency = Integer.parseInt(r);
                if (rescalingFrequency < 1) {
                    rescalingFrequency = RESCALE_FREQUENCY;
                }
            }

            String d = System.getProperty(DELAY_SCALING_PROPERTY);
            if (d != null) {
                this.delayRescalingUntilUnderflow = Boolean.parseBoolean(d);
            }

            if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
                if (stateCount == 4 && patternList.getPatternCount() < 10000) // TODO determine good cut-off
                    preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
            }

            boolean forceVectorization = false;
            String vectorizationString = System.getProperty(FORCE_VECTORIZATION);
            if (vectorizationString != null) {
                forceVectorization = true;
            }

            String tc = System.getProperty(THREAD_COUNT);
            if (tc != null) {
                threadCount = Integer.parseInt(tc);
            }

            if (threadCount < 1) {
                preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
                preferenceFlags |= BeagleFlag.THREADING_NONE.getMask();
            } else {
                preferenceFlags &= ~BeagleFlag.THREADING_NONE.getMask();
                preferenceFlags |= BeagleFlag.THREADING_CPP.getMask();
            }

            if (BeagleFlag.VECTOR_SSE.isSet(preferenceFlags) && (stateCount != 4)
                    && !forceVectorization && !IS_ODD_STATE_SSE_FIXED()
                    ) {
                // @todo SSE doesn't seem to work for larger state spaces so for now we override the
                // SSE option.
                preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
                preferenceFlags |= BeagleFlag.VECTOR_NONE.getMask();

                if (stateCount > 4 && this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                    this.rescalingScheme = PartialsRescalingScheme.DELAYED;
                }
            }

//            if (!BeagleFlag.PRECISION_SINGLE.isSet(preferenceFlags)) {
//                // if single precision not explicitly set then prefer double
//                preferenceFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
//            }

//            requirementFlags |= BeagleFlag.PRECISION_SINGLE.getMask();

            if (evolutionaryProcessDelegate.canReturnComplexDiagonalization()) {
                requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
            }

            if ((resourceList == null &&
                    (BeagleFlag.PROCESSOR_GPU.isSet(preferenceFlags) ||
                            BeagleFlag.FRAMEWORK_CUDA.isSet(preferenceFlags) ||
                            BeagleFlag.FRAMEWORK_OPENCL.isSet(preferenceFlags)))
                    ||
                    (resourceList != null && !IS_CPU_RESOURCE(resourceList[0]))) {
                // non-CPU implementations don't have SSE so remove default preference for SSE
                // when using non-CPU preferences or prioritising non-CPU resource
                preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
                preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
            }

            if (settings.usePreOrder && stateCount > 4) {
                requirementFlags |= BeagleFlag.PREORDER_TRANSPOSE_AUTO.getMask();
            }

            // after the default preferences above, which apply only while preferenceFlags is 0
            if (settings.useSpectralRepresentation && Boolean.getBoolean(SPECTRAL_CACHE_PROPERTY)) {
                preferenceFlags |= SPECTRAL_CACHE_FLAG;
            }
            if (Boolean.getBoolean(ADJOINT_BATCHED_PROPERTY)) {
                preferenceFlags |= ADJOINT_BATCHED_FLAG;
            }

            // start auto resource selection
            String resourceAuto = System.getProperty(RESOURCE_AUTO_PROPERTY);
            if (Boolean.parseBoolean(resourceAuto)) {

                long benchmarkFlags = 0;

                if (this.rescalingScheme == PartialsRescalingScheme.NONE) {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_NONE.getMask();
                } else if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS) {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_ALWAYS.getMask();
                } else {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_DYNAMIC.getMask();
                }

                logger.info("\t\tRunning benchmarks to automatically select fastest BEAGLE resource for analysis or partition... ");

                List<BenchmarkedResourceDetails> benchmarkedResourceDetails =
                        BeagleFactory.getBenchmarkedResourceDetails(
                                tipCount,
                                compactPartialsCount,
                                stateCount,
                                patternCount,
                                categoryCount,
                                resourceList,
                                preferenceFlags,
                                requirementFlags,
                                1, // eigenModelCount,
                                1, // partitionCount,
                                0, // calculateDerivatives,
                                benchmarkFlags);

                logger.info(" Benchmark results, from fastest to slowest:");

                for (BenchmarkedResourceDetails benchmarkedResource : benchmarkedResourceDetails) {
                    logger.info(benchmarkedResource.toString());
                }

                resourceList = new int[]{benchmarkedResourceDetails.get(0).getResourceNumber()};
            }
            // end auto resource selection


            beagle = BeagleFactory.loadBeagleInstance(
                    tipCount,
                    numPartials,
                    compactPartialsCount,
                    stateCount,
                    patternCount,
                    evolutionaryProcessDelegate.getEigenBufferCount(),
                    numMatrices,
                    categoryCount,
                    numScaleBuffers, // Always allocate; they may become necessary
                    resourceList,
                    preferenceFlags,
                    requirementFlags
            );

            if (growable) {
                try {
                    beagle.ensureBufferCounts(0, 0, 0); // every CPU instance from BEAGLE 4.2.0 can grow
                } catch (BeagleException e) {
                    if (beagle instanceof GeneralBeagleImpl) { // -java: it cannot do degree-2 operations either
                        throw new IllegalStateException("Degree-2 nodes at the epoch transition times need a " +
                                "native BEAGLE CPU instance, not the Java implementation (--java)", e);
                    }
                    throw new IllegalStateException("BEAGLE " + BeagleInfo.getVersion() + " cannot grow buffers (" +
                            e.getMessage() + "); check the native BEAGLE libraries or run with -D" +
                            PREALLOCATE_AUGMENTED_PROPERTY + "=true", e);
                }
            }

            InstanceDetails instanceDetails = beagle.getDetails();
            ResourceDetails resourceDetails = null;

            if (instanceDetails != null) {
                resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
                if (resourceDetails != null) {
                    StringBuilder sb = new StringBuilder("  Using BEAGLE resource ");
                    sb.append(resourceDetails.getNumber()).append(": ");
                    sb.append(resourceDetails.getName()).append("\n");
                    if (resourceDetails.getDescription() != null) {
                        String[] description = resourceDetails.getDescription().split("\\|");
                        for (String desc : description) {
                            if (desc.trim().length() > 0) {
                                sb.append("    ").append(desc.trim()).append("\n");
                            }
                        }
                    }
                    sb.append("    with instance flags: ").append(instanceDetails.toString());
                    logger.info(sb.toString());
                } else {
                    logger.info("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
                }
            } else {
                logger.info("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
            }

            instanceFlags = instanceDetails.getFlags();

            if ((instanceFlags & BeagleFlag.THREADING_CPP.getMask()) != 0) {
                if (IS_THREAD_COUNT_COMPATIBLE() && threadCount != 0) {
                    if (threadCount > 0) {
                        beagle.setCPUThreadCount(threadCount);
                        logger.info("    Using " + threadCount + " thread" + (threadCount > 1 ? "s" : "") + " for CPU.");
                    } else {
                        // if no thread_count is specified then this will be -1 so put no upper bound on threads
                        // currently the parser provides a default based on the number of cores as BEAGLE's
                        // default is suboptimal
                        logger.info("    Using default thread count for CPU.");
                        // this is just intended to remove the cap on number of threads so BEAGLE will
                        // make its own decision (for better or worse).
                        beagle.setCPUThreadCount(1000);
                    }
                } else {
                    logger.info("    BEAGLE threading turned off (or unavailable) for CPU.");
                }
            }

            if (patternList instanceof UncertainSiteList) { // TODO Remove
                useAmbiguities = true;
            }

//            if (!IS_PRE_ORDER_SUPPORTED() && settings.usePreOrder) {
//                throw new IllegalArgumentException("BEAGLE library does not support pre-order computation");
//            }

            //add in logger info for preOrder traversal
            logger.info("    " + (settings.usePreOrder ? "Using" : "Ignoring") + " preOrder partials in tree likelihood.");
            logger.info("    " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            logger.info("    With " + patternList.getPatternCount() + " unique site patterns.");

            if (patternList.areUncertain() && !useAmbiguities) {
                logger.info("    WARNING: Uncertain site patterns will be ignored.");
            }

            for (int i = 0; i < tipCount; i++) {
                // Find the id of tip i in the patternList
                String id = tree.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    if (isGhostTip == null || !isGhostTip[i]) {
                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + tree.getId() +
                                ", is not found in patternList, " + patternList.getId());
                    }
                    // all unknown, so that it can stand in for a ghost subtree when one-child operations are not used
                    setGhostTip(beagle, i, useAmbiguities);
                } else {
                    if (useAmbiguities) {
                        setPartials(beagle, patternList, index, i);
                    } else {
                        setStates(beagle, patternList, index, i);
                    }
                }
            }

            beagle.setPatternWeights(patternWeights);

            if (isGhostTip != null) {
                setUpGhosts(logger);
            }

            if (epochDelegate != null) {
                logger.info("    Adding degree-2 nodes at the epoch transition times (up to "
                        + augmentedNodeCapacity + " extra nodes, BEAGLE buffers "
                        + (growable ? "grown on demand" : "preallocated") + ").");
            }

            String rescaleMessage = "    Using rescaling scheme : " + this.rescalingScheme.getText();
            if (this.rescalingScheme == PartialsRescalingScheme.AUTO &&
                    resourceDetails != null &&
                    (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
                // If auto scaling in BEAGLE is not supported then do it here
                this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                rescaleMessage = "    Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText();
            }
            boolean parenthesis = false;
            if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                rescaleMessage += " (rescaling every " + rescalingFrequency + " evaluations";
                parenthesis = true;
            }
            if (this.delayRescalingUntilUnderflow) {
                rescaleMessage += (parenthesis ? ", " : "(") + "delay rescaling until first overflow";
                parenthesis = true;
            }
            rescaleMessage += (parenthesis ? ")" : "");
            logger.info(rescaleMessage);

            if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                everUnderflowed = false; // If false, BEAST does not rescale until first under-/over-flow.
            }

            updateSubstitutionModel = true;
            updateSiteModel = true;
            updateRootFrequency = true;

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }

        instanceCount++;
    }

    public BeagleDataLikelihoodDelegate(Tree tree,
                                        PatternList patternList,
                                        BranchModel branchModel,
                                        SiteRateModel siteRateModel,
                                        boolean useAmbiguities,
                                        boolean preferGPU,
                                        PartialsRescalingScheme rescalingScheme,
                                        boolean delayRescalingUntilUnderflow,
                                        PreOrderSettings settings

    ) {
        this(tree, patternList, branchModel, siteRateModel, useAmbiguities, preferGPU, rescalingScheme,
                delayRescalingUntilUnderflow, settings, null, null);
    }

    @Override
    public String getReport() {
        return null;
    }

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        if ((instanceFlags & BeagleFlag.FRAMEWORK_CPU.getMask()) != 0) {
            return TreeTraversal.TraversalType.POST_ORDER;
        } else {
            return TreeTraversal.TraversalType.REVERSE_LEVEL_ORDER;
        }
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
    public RateRescalingScheme getRateRescalingScheme() {
        return RateRescalingScheme.NONE;
    }

    public final BranchModel getBranchModel() {
        return branchModel;
    }

    public PatternList getPatternList() {
        return this.patternList;
    }

    public Beagle getBeagleInstance() { return beagle; }

    private int getSingleScaleBufferCount() {
        return internalNodeCount + 1;
    }

    // BEAGLE buffer counts with buffers for a augmented nodes; only with augmented epoch nodes

    private int partialsFor(int a) {
        // BEAGLE's index space is the partials count plus the compact count, and the helper's count already
        // includes the tips, so without this subtraction tipCount internal buffers would never be used
        return augmentedPartialHelper.getBufferCount(a) - compactPartialsCount;
    }

    private int matricesFor(int a) {
        return epochDelegate.getMatrixBufferCount(a);
    }

    private int scaleFor(int a) {
        // with pre-order, one gradient scale buffer per original and augmented node, after the others
        return scaleBufferHelper.getBufferCount() + (settings.usePreOrder ? nodeCount + a : 0);
    }

    private void setAugmentedAllocation(int a) {
        augmentedPartialHelper.setAllocated(a);
        epochDelegate.setAugmentedBufferCount(a);
        allocatedAugmented = a;
    }

    /**
     * Gives BEAGLE buffers to every augmented node the registry has used, before any BEAGLE call of the evaluation:
     * every id is assigned by the traversal that precedes calculateLikelihood. Exactly as many as needed the first
     * time, then at least 1.5 times as many, at most the capacity of the registry. Indices never move.
     */
    private void ensureAugmentedBuffers() {
        final int needed = registry.getHighWaterMark();
        if (needed <= allocatedAugmented) {
            return;
        }
        if (!growable) { // unreachable: the registry never exceeds its capacity, which is preallocated
            throw new IllegalStateException("Augmented epoch nodes exceed their BEAGLE buffers");
        }
        final int previous = allocatedAugmented;
        final int target = Math.min(registry.getCapacity(), Math.max(needed, previous + (previous >> 1)));
        beagle.ensureBufferCounts(partialsFor(target), matricesFor(target), scaleFor(target)); // throws on failure
        setAugmentedAllocation(target);
        ++growthEventCount;
        Logger.getLogger("dr.evomodel").fine("Degree-2 buffer growth" + (getId() != null ? " (" + getId() + ")" : "") +
                ": " + previous + " -> " + target + " augmented nodes, " + needed + " used");
    }

    /**
     * @return the number of times the BEAGLE buffers of augmented nodes grew
     */
    public int getGrowthEventCount() {
        return growthEventCount;
    }

    /**
     * @return the number of augmented nodes with BEAGLE buffers
     */
    public int getAllocatedAugmentedCount() {
        return allocatedAugmented;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternList   patternList
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    private final void setPartials(Beagle beagle,
                                   PatternList patternList,
                                   int sequenceIndex,
                                   int nodeIndex) {
        double[] partials = new double[patternCount * stateCount * categoryCount];

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

            if (patternList instanceof UncertainSiteList) {
                ((UncertainSiteList) patternList).fillPartials(sequenceIndex, i, partials, v);
                v += stateCount;
                // TODO Add this functionality to SimpleSiteList to avoid if statement here
            } else if (patternList.areUncertain()) {

                double[] prob = patternList.getUncertainPatternState(sequenceIndex, i);
                System.arraycopy(prob, 0, partials, v, stateCount);
                v += stateCount;

            } else {
                int state = patternList.getPatternState(sequenceIndex, i);
                boolean[] stateSet = dataType.getStateSet(state);

                for (int j = 0; j < stateCount; j++) {
                    if (stateSet[j]) {
                        partials[v] = 1.0;
                    } else {
                        partials[v] = 0.0;
                    }
                    v++;
                }
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * stateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(nodeIndex, partials);
    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    private final void setPartials(Beagle beagle,
                                   TipStatesModel tipStatesModel,
                                   int nodeIndex) {
        double[] partials = new double[patternCount * stateCount * categoryCount];

        tipStatesModel.getTipPartials(nodeIndex, partials);

        // if there is more than one category then replicate the partials for each
        int n = patternCount * stateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(nodeIndex, partials);
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternList   patternList
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    private final void setStates(Beagle beagle,
                                 PatternList patternList,
                                 int sequenceIndex,
                                 int nodeIndex) {
        int i;

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {

            states[i] = patternList.getPatternState(sequenceIndex, i);
        }

        beagle.setTipStates(nodeIndex, states);
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) throws LikelihoodException {

        if (registry != null) {
            ensureAugmentedBuffers();
        }

        if(siteAssignInd != null) {
            setAllPatternWeights(siteAssignInd);
            beagle.setPatternWeights(patternWeights);
        }

        //recomputeScaleFactors = false;
        if (DEBUG) {
            System.out.println("Partition: " + this.getModelName());
        }

        if (!this.delayRescalingUntilUnderflow || everUnderflowed) {
            if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS || this.rescalingScheme == PartialsRescalingScheme.DELAYED) {
                useScaleFactors = true;
                recomputeScaleFactors = true;
            } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                useScaleFactors = true;

                if (rescalingCount > rescalingFrequency) {
                    if (DEBUG) {
                        System.out.println("rescalingCount > rescalingFrequency");
                    }
                    rescalingCount = 0;
                    rescalingCountInner = 0;
                }

                if (DEBUG) {
                    System.out.println("rescalingCountInner = " + rescalingCountInner);
                }

                if (rescalingCountInner < RESCALE_TIMES) {
                    if (DEBUG) {
                        System.out.println("rescalingCountInner < RESCALE_TIMES");
                    }

                    recomputeScaleFactors = true;

                    rescalingCountInner++;

                    throw new LikelihoodRescalingException();

                }

                //underflowHandling takes into account the first evaluation when initiating the MCMC chain
                //suggest replacing with boolean initialEvaluation
                if (initialEvaluation) {
                    if (underflowHandling < 1) {
                        underflowHandling++;
                        if (DEBUG) {
                            System.out.println("underflowHandling < 1");
                        }
                    } else if (underflowHandling == 1) {
                        if (DEBUG) {
                            System.out.println("underflowHandling == 1");
                        }
                        recomputeScaleFactors = true;
                        underflowHandling++;
                        initialEvaluation = false;
                    }
                }

                rescalingCount++;
            }
        }

        if (RESCALING_OFF) { // a debugging switch
            useScaleFactors = false;
            recomputeScaleFactors = false;
        }

        int becameRealCount = 0;
        if (nodeStatus != null) {
            evaluationStamp++;
            rootInputsChanged = false;
            becameRealCount = updateGhostStatus(nodeOperations, rootNodeNumber);
        }

        int branchUpdateCount = 0;
        if (nodeStatus == null) {
            for (BranchOperation op : branchOperations) {
                branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
                branchLengths[branchUpdateCount] = op.getBranchLength();
                branchUpdateCount++;
            }
        } else if (mergeGhostBranches) {
            branchUpdateCount = collectMergedBranches(branchOperations, nodeOperations);
        } else {
            branchUpdateCount = collectRealBranches(branchOperations, becameRealCount, rootNodeNumber);
        }

        if (updateSubstitutionModel) { // TODO More efficient to update only the substitution model that changed, instead of all
            evolutionaryProcessDelegate.updateSubstitutionModels(beagle, flip);

            // we are currently assuming a no-category model...
        }

        if (updateSiteModel) {
            double[] categoryRates = this.siteRateModel.getCategoryRates();
            if (categoryRates == null) {
                // If this returns null then there was a numerical error calculating the category rates
                // (probably a very small alpha) so reject the move.
                return Double.NEGATIVE_INFINITY;
            }
            beagle.setCategoryRates(categoryRates);
            // TODO Try beagle.setCategoryWeights() here
            double[] categoryWeights = this.siteRateModel.getCategoryProportions();
            if (mergeGhostBranches) {
                rootInputsChanged |= !Arrays.equals(categoryRates, lastCategoryRates) ||
                        !Arrays.equals(categoryWeights, lastCategoryWeights);
                lastCategoryRates = categoryRates.clone();
                lastCategoryWeights = categoryWeights.clone();
            }
            // these could be set only when they change but store/restore would need to be considered
            beagle.setCategoryWeights(0, categoryWeights); // TODO move
        }

        if (updateRootFrequency) {
            double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
            beagle.setStateFrequencies(0, frequencies); // TODO make lazy?
            if (mergeGhostBranches) {
                rootInputsChanged |= !Arrays.equals(frequencies, lastRootFrequencies);
                lastRootFrequencies = frequencies.clone();
            }
        }

        if (branchUpdateCount > 0) {
            evolutionaryProcessDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount,
                    flip);
        }

        if (COUNT_CALCULATIONS) {
            totalMatrixUpdateCount += branchUpdateCount;
        }

        if (flip) {
            // Flip all the buffers to be written to first...
            for (NodeOperation op : nodeOperations) {
                final int nodeNum = op.getNodeNumber();
                if (nodeStatus == null || (mergeGhostBranches ? isRecomputed(nodeNum) : isComputed(nodeStatus[nodeNum]))) {
                    partialBufferHelper.flipOffset(nodeNum);
                }
            }
        }

        int operationCount = 0;
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            int nodeNum = op.getNodeNumber();

            int leftChild = op.getLeftChild();
            int rightChild = op.getRightChild();
            int rightMatrixNode = rightChild;
            boolean fromOneChild = false;

            if (nodeStatus != null) {
                final byte status = nodeStatus[nodeNum];
                if (status == GHOST || status == PASS) {
                    // nothing below a ghost has any data, so it contributes nothing and is not computed; a node
                    // passed through is not computed either, its parent reading its real child across both branches
                    continue;
                }
                if (mergeGhostBranches && !isRecomputed(nodeNum)) {
                    // none of its inputs changed, so it keeps its partials, and its scale factor with them
                    continue;
                }
                if (mergeGhostBranches) {
                    // the nearest descendants with data, each across its joined branch
                    leftChild = effectiveLeft[nodeNum];
                    rightChild = effectiveRight[nodeNum];
                    rightMatrixNode = rightChild;
                }
                if (status == ONE_CHILD) {
                    final int realChild = mergeGhostBranches ? leftChild :
                            ((nodeStatus[leftChild] == GHOST) ? rightChild : leftChild);
                    leftChild = realChild;
                    if (useOneChildOperations) {
                        // computed from the real child alone, but rescaled like any other node: a lineage can pass
                        // through a long run of such nodes, and the partials would drift towards underflow
                        fromOneChild = true;
                    } else {
                        // an all-unknown tip stands in for the ghost subtree: its contribution is one whatever the
                        // matrix, so the real child's matrix is used for it too
                        rightChild = ghostTipNumber;
                        rightMatrixNode = realChild;
                    }
                }
            }

            operationCount++;

            operations[k] = partialBufferHelper.getOffsetIndex(nodeNum);

            if (rightChild == NodeOperation.NO_CHILD && !fromOneChild) {
                // A degree-2 node has one child and is never rescaled
                setDegree2Operation(k, leftChild);
                k += Beagle.OPERATION_TUPLE_SIZE;
                continue;
            }

            if (!isRestored && !partialBufferHelper.isSafeUpdate(nodeNum) && !recomputeScaleFactors) {
                System.err.println("Stored partial should not be updated!");
            }

            if (useScaleFactors) {
                // get the index of this scaling buffer
                int n = nodeNum - tipCount;

                if (recomputeScaleFactors) {
                    // flip the indicator: can take either n or (internalNodeCount + 1) - n
                    scaleBufferHelper.flipOffset(n);

                    // store the index
                    scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                    operations[k + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                    operations[k + 2] = Beagle.NONE;

                    if (nodeStatus != null) {
                        hasScaleFactor[nodeNum] = true;
                    }

                } else if (nodeStatus != null && !hasScaleFactor[nodeNum]) {
                    // a node that has started to be computed since the factors were last recomputed has none, and
                    // is computed without one until they next are, as writing one now would fit it to partials
                    // that have drifted, and that drift would be handed on if it later lost it
                    operations[k + 1] = Beagle.NONE;
                    operations[k + 2] = Beagle.NONE;

                } else {
                    operations[k + 1] = Beagle.NONE;
                    operations[k + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                }

            } else {

                if (useAutoScaling) {
                    scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                }
                operations[k + 1] = Beagle.NONE; // Not using scaleFactors
                operations[k + 2] = Beagle.NONE;
            }

            operations[k + 3] = partialBufferHelper.getOffsetIndex(leftChild); // source node 1
            operations[k + 4] = evolutionaryProcessDelegate.getMatrixIndex(leftChild); // source matrix 1
            if (fromOneChild) {
                operations[k + 5] = Beagle.NONE; // no source node 2
                operations[k + 6] = Beagle.NONE;
            } else {
                operations[k + 5] = partialBufferHelper.getOffsetIndex(rightChild); // source node 2
                operations[k + 6] = evolutionaryProcessDelegate.getMatrixIndex(rightMatrixNode); // source matrix 2
            }

            k += Beagle.OPERATION_TUPLE_SIZE;
        }

        if (mergeGhostBranches && operationCount == 0 && !isFullRescale() && leavingScaledCount == 0 &&
                !rootInputsChanged && rootNodeNumber == lastRootNode && !Double.isNaN(lastLogL) && siteAssignInd == null) {
            // Nothing that is computed changed, so neither would the likelihood: a move that touched ghosts alone.
            updateSubstitutionModel = false;
            updateSiteModel = false;
            updateRootFrequency = false;
            return lastLogL;
        }

        final boolean patchScaleTotal = nodeStatus != null && useScaleFactors && !recomputeScaleFactors &&
                leavingScaledCount > 0;
        if (patchScaleTotal) {
            removeOldScaleFactors();
        }

        if (useScaleFactors && recomputeScaleFactors) {
            ++scaleFactorsVersion; // this call writes individual scale factors (see isNegAncestorScaleBufferCurrent)
        }
        beagle.updatePartials(operations, operationCount, Beagle.NONE);

        if (COUNT_CALCULATIONS) {
            totalEvaluationCount += 1;
            totalPartialsUpdateCount += operationCount;
        }

        int rootIndex = partialBufferHelper.getOffsetIndex(rootNodeNumber);

        int cumulateScaleBufferIndex = Beagle.NONE;

        int[] accumulatedScaleBuffers = scaleBufferIndices;
        int accumulatedScaleBufferCount = internalNodeCount;
        if (nodeStatus != null && ((useScaleFactors && recomputeScaleFactors) || useAutoScaling)) {
            accumulatedScaleBufferCount = collectScaledNodes();
            accumulatedScaleBuffers = scaledNodeBufferIndices;
        }

        if (useScaleFactors) {
            if (recomputeScaleFactors) {
                scaleBufferHelper.flipOffset(internalNodeCount);
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                beagle.resetScaleFactors(cumulateScaleBufferIndex);
                beagle.accumulateScaleFactors(accumulatedScaleBuffers, accumulatedScaleBufferCount,
                        cumulateScaleBufferIndex);
            } else {
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
            }
        } else if (useAutoScaling) {
            beagle.accumulateScaleFactors(accumulatedScaleBuffers, accumulatedScaleBufferCount, Beagle.NONE);
        }

        double[] sumLogLikelihoods = new double[1];

        if (DEBUG) {
            System.out.println("useScaleFactors=" + useScaleFactors + " recomputeScaleFactors=" + recomputeScaleFactors + " (" + getId() + ")");
        }

        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

        double logL = sumLogLikelihoods[0];

        /*if (DEBUG) {
            System.out.println(logL);
            if (logL > -90000) {
                System.exit(0);
            }
        }*/

        if (Double.isNaN(logL) || Double.isInfinite(logL)) {

            if (DEBUG) {
                System.out.println("Double.isNaN(logL) || Double.isInfinite(logL) (" + getId() + ")");
            }

            everUnderflowed = true;

            logL = Double.NEGATIVE_INFINITY;

            if (firstRescaleAttempt && (delayRescalingUntilUnderflow || rescalingScheme == PartialsRescalingScheme.DELAYED)) {

                if (rescalingScheme == PartialsRescalingScheme.DYNAMIC || (rescalingCount == 0)) {
                    // show a message but only every 1000 rescales
                    if (rescalingMessageCount % 1000 == 0) {
                        if (rescalingMessageCount > 0) {
                            Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood (" + rescalingMessageCount + " messages not shown; " + getId() + ").");
                        } else {
                            if (getId() != null) {
                                Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling... (" + getId() + ")");
                            } else {
                                Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling...");
                            }
                        }
                    }
                    rescalingMessageCount += 1;
                }

                useScaleFactors = true;
                recomputeScaleFactors = true;

                firstRescaleAttempt = false; // Only try to rescale once

                rescalingCount--;

            }

            // Double buffering stays on for the retry. A buffer is flipped only once before an accept or reject, so
            // the buffers this attempt wrote are simply written again, while every other node, which the retry also
            // recomputes, writes to its spare buffer. Turning flipping off instead overwrote the stored partials of
            // those nodes with partials rescaled by new factors: if the proposal was then rejected, the old scale
            // factors were restored over the new partials and the likelihood of the current state was wrong.
            underflowHandling = 0;

            if (USE_CACHED_EXCEPTION) {
                if (cachedException == null) {
                    cachedException = new LikelihoodUnderflowException();
                }
                throw cachedException;
            } else {
                throw new LikelihoodUnderflowException();
            }

        } else {

            firstRescaleAttempt = true;
            recomputeScaleFactors = false;
            flip = true;
            if (mergeGhostBranches) {
                lastLogL = logL;
                lastRootNode = rootNodeNumber;
            }

        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
        updateRootFrequency = false;
        //********************************************************************

        // If these are needed...
        //if (patternLogLikelihoods == null) {
        //    patternLogLikelihoods = new double[patternCount];
        //}
        //beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        return logL;
    }

    public double[] getSiteLogLikelihoods(){
        double[] patternLogLikelihoods = new double[patternCount];
        beagle.getSiteLogLikelihoods(patternLogLikelihoods);
        return patternLogLikelihoods;
    }

    public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        /* No need to rescale partials */
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
    }

    private void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }

    private void setAllPatternWeights(Parameter siteAssignInd){
        // First, set all patternWeights to 0.0
        Arrays.fill(patternWeights, 0.0);
        for(int i = 0; i < siteAssignInd.getSize(); i++){
            int cat = (int) siteAssignInd.getParameterValue(i);
            if((int) partitionCat.getParameterValue(0) == cat) {
                int pwIndex = patternList.getPatternIndex(i);
                patternWeights[pwIndex] = patternWeights[pwIndex] + 1;
            }
        }
    }

    public int getPartitionCat(){
        return (int) partitionCat.getParameterValue(0);
    };

    @Override
    public void makeDirty() {
        updateSiteModel = true;
        updateSubstitutionModel = true;
        updateRootFrequency = true;
        fireModelChanged();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == siteRateModel) {
            updateSiteModel = true;
        } else if (model == branchModel || (epochTimeProvider != null && model != epochTimeProvider)) {
            // with augmented epoch nodes, the substitution models and root frequencies stand in for the branch model
            updateSubstitutionModel = true;
            updateRootFrequency = true;
        }

        // Tell TreeDataLikelihood to update all nodes
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable == siteAssignInd) {
            setAllPatternWeights(siteAssignInd);
            beagle.setPatternWeights(patternWeights);
        }
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void storeState() {
        if(siteAssignInd != null) {
            System.arraycopy(patternWeights, 0, storedPatternWeights, 0, patternWeights.length);
        }

        partialBufferHelper.storeState();
        evolutionaryProcessDelegate.storeState();

        if (nodeStatus != null) {
            System.arraycopy(nodeStatus, 0, storedNodeStatus, 0, nodeStatus.length);
            System.arraycopy(hasScaleFactor, 0, storedHasScaleFactor, 0, hasScaleFactor.length);
            System.arraycopy(ghostBranchLengths, 0, storedGhostBranchLengths, 0, ghostBranchLengths.length);
            if (mergeGhostBranches) {
                System.arraycopy(passChild, 0, storedPassChild, 0, passChild.length);
                System.arraycopy(mergedLength, 0, storedMergedLength, 0, mergedLength.length);
                System.arraycopy(effectiveLeft, 0, storedEffectiveLeft, 0, effectiveLeft.length);
                System.arraycopy(effectiveRight, 0, storedEffectiveRight, 0, effectiveRight.length);
                storedPreviousRoot = previousRoot;
                storedLastLogL = lastLogL;
                storedLastRootNode = lastRootNode;
                storedLastCategoryRates = lastCategoryRates;      // replaced, never changed in place
                storedLastCategoryWeights = lastCategoryWeights;
                storedLastRootFrequencies = lastRootFrequencies;
            }
        }

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }

        // turn on double buffering flipping (may have been turned off to enable a rescale)
        flip = true;

        isRestored = false;
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restoreState() {

        if (siteAssignInd != null) {
            //updateSiteModel = true;
            updateSubstitutionModel = true;
            // Swap pointers
            double[] temp = patternWeights;
            patternWeights = storedPatternWeights;
            storedPatternWeights = temp;
        }

        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore
        updateRootFrequency = true;

        partialBufferHelper.restoreState();
        evolutionaryProcessDelegate.restoreState();

        if (nodeStatus != null) {
            byte[] tmpStatus = storedNodeStatus;
            storedNodeStatus = nodeStatus;
            nodeStatus = tmpStatus;

            boolean[] tmpScaled = storedHasScaleFactor;
            storedHasScaleFactor = hasScaleFactor;
            hasScaleFactor = tmpScaled;

            double[] tmpLengths = storedGhostBranchLengths;
            storedGhostBranchLengths = ghostBranchLengths;
            ghostBranchLengths = tmpLengths;

            if (mergeGhostBranches) {
                int[] tmpChild = storedPassChild;
                storedPassChild = passChild;
                passChild = tmpChild;

                double[] tmpMerged = storedMergedLength;
                storedMergedLength = mergedLength;
                mergedLength = tmpMerged;

                int[] tmpEffective = storedEffectiveLeft;
                storedEffectiveLeft = effectiveLeft;
                effectiveLeft = tmpEffective;
                tmpEffective = storedEffectiveRight;
                storedEffectiveRight = effectiveRight;
                effectiveRight = tmpEffective;

                previousRoot = storedPreviousRoot;
                lastLogL = storedLastLogL;
                lastRootNode = storedLastRootNode;
                lastCategoryRates = storedLastCategoryRates;
                lastCategoryWeights = storedLastCategoryWeights;
                lastRootFrequencies = storedLastRootFrequencies;
            }
        }

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp;
//            rescalingCount = storedRescalingCount;
        }

        isRestored = true;

    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
        // Do nothing
    }

    /**
     * Fills in the sources of an operation for a node with a single child (the destination is already set); the
     * second child is Beagle.NONE. Rescaling the children keeps the scale factors valid because propagating along a
     * segment by a stochastic matrix cannot make partials any smaller.
     */
    private void setDegree2Operation(int k, int child) {
        operations[k + 1] = Beagle.NONE; // Not writing scaleFactors
        operations[k + 2] = Beagle.NONE; // Not reading scaleFactors

        operations[k + 3] = partialBufferHelper.getOffsetIndex(child); // source node 1
        operations[k + 4] = evolutionaryProcessDelegate.getMatrixIndex(child); // source matrix 1

        operations[k + 5] = Beagle.NONE; // no source node 2
        operations[k + 6] = Beagle.NONE;
    }

    @Override
    public void setComputePostOrderStatisticsOnly(boolean computePostOrderStatistic) {
        // Do nothing
    }

    @Override
    public boolean providesPostOrderStatisticsOnly() { return false; }

    @Override
    public int vectorizeNodeOperations(List<ProcessOnTreeDelegate.NodeOperation> nodeOperations, int[] operations) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void acceptState() {
    }

    public final boolean isUsePreOrder(){
        return this.settings.usePreOrder;
    }

    public final EvolutionaryProcessDelegate getEvolutionaryProcessDelegate() {
        return this.evolutionaryProcessDelegate;
    }

    public final SiteRateModel getSiteRateModel(){
        return this.siteRateModel;
    }

    public final int getPartialBufferIndex(int nodeNumber) {
        return partialBufferHelper.getOffsetIndex(nodeNumber);
    }

    public final int getScaleBufferCount() {
        return scaleBufferHelper.getBufferCount();
    }

    /**
     * @return whether the last likelihood evaluation actually applied rescaling
     * (independent of the requested PartialsRescalingScheme -- e.g. underflow
     * handling or a scheme-parsing fallback can still trigger it).
     */
    public final boolean isUsingScaleFactors() {
        return useScaleFactors;
    }

    /**
     * Asks the next likelihood evaluation for fresh scale factors at the nodes it updates. The dynamic scheme
     * otherwise reuses those of an earlier evaluation; after a large change (a different tree, say), partials scaled
     * with them can grow so large that a pre-order traversal from them overflows, though the likelihood is finite.
     *
     * @return false if the scale factors are always fresh or there are none, so that asking changes nothing
     */
    public final boolean recomputeScaleFactorsNextEvaluation() {
        if (rescalingScheme != PartialsRescalingScheme.DYNAMIC || !useScaleFactors) {
            return false;
        }
        recomputeScaleFactors = true;
        return true;
    }

    /**
     * @return the scale buffer holding internal node nodeNumber's own individual
     * (not cumulative) scale factor from the last likelihood evaluation, or
     * Beagle.NONE if nodeNumber is a tip or rescaling was not applied. Callers
     * needing a node's full CUMULATIVE post-order scale (its own factor plus
     * every rescaled descendant's) must build that themselves, e.g. via
     * beagleAccumulateScaleFactors over every rescaled node at-or-below it --
     * see beagleCalculateAdjointCrossProductDerivative's contract in beagle.h.
     */
    public final int getNodeIndividualScaleBufferIndex(int nodeNumber) {
        if (!useScaleFactors || nodeNumber < tipCount || nodeNumber >= nodeCount) { // not augmented nodes
            return Beagle.NONE;
        }
        return scaleBufferIndices[nodeNumber - tipCount];
    }

    /**
     * Whether gradient scale buffer slot (from getScaleBufferCount() on) still holds what
     * AbstractBeagleGradientDelegate.buildNegAncestorCumulativeScaleBuffer would build into it from source (a slot, or
     * Beagle.NONE for a reset) less individual (a node's individual scale buffer): true if it was last built from the
     * same two, source has not been rebuilt since, and no individual scale factors have been written since. Every
     * delegate on this likelihood shares the slots, so this record of them lives here, with the BEAGLE instance.
     */
    public final boolean isNegAncestorScaleBufferCurrent(int slot, int source, int individual) {
        final int i = slot - getScaleBufferCount();
        if (i < 0 || i >= negAncestorBuild.length || negAncestorBuild[i] == 0) {
            return false;
        }
        final long sourceBuild;
        if (source == Beagle.NONE) {
            sourceBuild = 0;
        } else {
            final int j = source - getScaleBufferCount();
            if (j < 0 || j >= negAncestorBuild.length) {
                return false;
            }
            sourceBuild = negAncestorBuild[j];
        }
        return negAncestorSource[i] == source && negAncestorIndividual[i] == individual &&
                negAncestorVersion[i] == scaleFactorsVersion && negAncestorSourceBuild[i] == sourceBuild;
    }

    /**
     * Records that slot was just built from source less individual (see isNegAncestorScaleBufferCurrent).
     */
    public final void markNegAncestorScaleBufferBuilt(int slot, int source, int individual) {
        final int i = slot - getScaleBufferCount();
        if (i < 0) {
            return;
        }
        if (i >= negAncestorBuild.length) {
            final int n = Math.max(i + 1, 2 * negAncestorBuild.length);
            negAncestorSource = Arrays.copyOf(negAncestorSource, n);
            negAncestorIndividual = Arrays.copyOf(negAncestorIndividual, n);
            negAncestorVersion = Arrays.copyOf(negAncestorVersion, n);
            negAncestorSourceBuild = Arrays.copyOf(negAncestorSourceBuild, n);
            negAncestorBuild = Arrays.copyOf(negAncestorBuild, n);
        }
        final int j = source - getScaleBufferCount();
        negAncestorSource[i] = source;
        negAncestorIndividual[i] = individual;
        negAncestorVersion[i] = scaleFactorsVersion;
        negAncestorSourceBuild[i] = (source == Beagle.NONE || j < 0 || j >= negAncestorBuild.length) ? 0 : negAncestorBuild[j];
        negAncestorBuild[i] = ++negAncestorBuildCount;
    }


    public final int getPartialBufferCount() {
        return partialBufferHelper.getBufferCount();
    }

    /**
     * @param id an original or augmented node
     * @return the buffer of its pre-order partial
     */
    public final int getPreOrderPartialIndex(int id) {
        if (augmentedPartialHelper == null) {
            return getPartialBufferCount() + id;
        }
        if (!settings.usePreOrder) { // these indices would be the post-order partials of augmented nodes
            throw new IllegalStateException("Pre-order partials with augmented epoch nodes need usePreOrder=\"true\"");
        }
        return (id < nodeCount) ? (2 * nodeCount - tipCount) + id : augmentedPartialHelper.getSlotIndex(id, 2);
    }

    public PreOrderSettings getPreOrderSettings() {
        return settings;
    }


    public boolean getPreferGPU() {
        return prefersGPU;
    }


    public boolean getUseAmbiguities() {
        return useAmbiguities;
    }

    public PartialsRescalingScheme getRescalingScheme() {
        return rescalingScheme;
    }

    public boolean getDelayRescalingUntilUnderflow() {
        return delayRescalingUntilUnderflow;
    }
    // **************************************************************
    // INSTANCE PROFILEABLE
    // **************************************************************

    @Override
    public long getTotalCalculationCount() {
        // Can only return one count at the moment so return the number of partials updated
        return totalPartialsUpdateCount;
    }

    // **************************************************************
    // INSTANCE CITABLE
    // **************************************************************

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "BEAGLE likelihood calculation library";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.AYRES_2019_BEAGLE);
    }

    private void releaseBeagle() throws Throwable {
        if (beagle != null && releaseSingleton) {
            beagle.finalize();
            releaseSingleton = false;
        }
    }

    public static void releaseBeagleDataLikelihoodDelegate(TreeDataLikelihood treeDataLikelihood) throws Throwable {
        DataLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (likelihoodDelegate instanceof BeagleDataLikelihoodDelegate) {
            BeagleDataLikelihoodDelegate delegate = (BeagleDataLikelihoodDelegate) likelihoodDelegate;
            delegate.releaseBeagle();
        }
    }

    public static void releaseAllBeagleInstances() throws Throwable {
        // every chain of a multi-chain run has its own registry
        for (dr.inference.model.ModelRegistry registry : dr.inference.model.ModelRegistry.getAll()) {
            for (Likelihood likelihood : registry.getFullLikelihoodSet()) {
                if (likelihood instanceof TreeDataLikelihood) {
                    releaseBeagleDataLikelihoodDelegate((TreeDataLikelihood) likelihood);
                } else if (likelihood instanceof CompoundLikelihood) {
                    for (Likelihood likelihood2: ((CompoundLikelihood) likelihood).getLikelihoods()) {
                        if (likelihood2 instanceof TreeDataLikelihood) {
                            releaseBeagleDataLikelihoodDelegate((TreeDataLikelihood) likelihood2);
                        }
                    }
                }
            }
        }
    }

    // **************************************************************
    // GHOST TAXA
    // **************************************************************

    /**
     * @return for each tip of the tree, whether its taxon is missing from the pattern list, or null if none is
     */
    private static boolean[] findGhostTips(Tree tree, PatternList patternList) {
        final int tipCount = tree.getExternalNodeCount();
        boolean[] ghosts = null;
        for (int i = 0; i < tipCount; i++) {
            if (patternList.getTaxonIndex(tree.getTaxonId(i)) == -1) {
                if (ghosts == null) {
                    ghosts = new boolean[tipCount];
                }
                ghosts[i] = true;
            }
        }
        return ghosts;
    }

    private static int countGhosts(boolean[] isGhostTip) {
        int count = 0;
        for (boolean ghost : isGhostTip) {
            if (ghost) {
                count++;
            }
        }
        return count;
    }

    /**
     * Gives a ghost tip data that is unknown at every site: a state outside the state space, or partials of one.
     */
    private void setGhostTip(Beagle beagle, int nodeIndex, boolean useAmbiguities) {
        if (useAmbiguities) {
            double[] partials = new double[patternCount * stateCount * categoryCount];
            Arrays.fill(partials, 1.0);
            beagle.setPartials(nodeIndex, partials);
        } else {
            int[] states = new int[patternCount];
            Arrays.fill(states, dataType.getUnknownState());
            beagle.setTipStates(nodeIndex, states);
        }
    }

    /**
     * Sets up the status of every node and decides how a node with one ghost child is computed, now that the BEAGLE
     * instance is known.
     */
    private void setUpGhosts(Logger logger) {
        nodeStatus = new byte[nodeCount];       // internal nodes start as FULL, and are set by the first evaluation
        storedNodeStatus = new byte[nodeCount];
        ghostBranchLengths = new double[nodeCount];
        storedGhostBranchLengths = new double[nodeCount];
        becameReal = new int[nodeCount];
        branchListed = new boolean[nodeCount];
        scaledNodeBufferIndices = new int[internalNodeCount];
        leavingScaled = new int[nodeCount];
        hasScaleFactor = new boolean[nodeCount];
        storedHasScaleFactor = new boolean[nodeCount];

        for (int i = 0; i < tipCount; i++) {
            if (isGhostTip[i]) {
                nodeStatus[i] = GHOST;
                if (ghostTipNumber < 0) {
                    ghostTipNumber = i;
                }
            }
        }
        System.arraycopy(nodeStatus, 0, storedNodeStatus, 0, nodeCount);

        // Joining the branches through a node with one ghost child needs one rate matrix for every branch, so that
        // a branch's matrix depends on its length alone.
        String mergeProblem = ghostMergeProblem(branchModel, settings);
        if (mergeProblem == null && !(evolutionaryProcessDelegate instanceof HomogenousSubstitutionModelDelegate)) {
            mergeProblem = "the transition matrices are not made from a single substitution model";
        }
        mergeGhostBranches = mergeProblem == null;
        if (mergeGhostBranches) {
            passChild = new int[nodeCount];
            storedPassChild = new int[nodeCount];
            mergedLength = new double[nodeCount];
            storedMergedLength = new double[nodeCount];
            Arrays.fill(mergedLength, Double.NaN);
            Arrays.fill(storedMergedLength, Double.NaN);
            effectiveLeft = new int[nodeCount];
            storedEffectiveLeft = new int[nodeCount];
            storedEffectiveRight = new int[nodeCount];
            recomputedAt = new int[nodeCount];
            statusChangedAt = new int[nodeCount];
            effectiveRight = new int[nodeCount];
            Arrays.fill(effectiveLeft, -2);
            Arrays.fill(effectiveRight, -2);
            Arrays.fill(storedEffectiveLeft, -2);
            Arrays.fill(storedEffectiveRight, -2);
        }

        // One-child operations (a second child of Beagle.NONE) are only in BEAGLE's CPU implementations from 4.2.0,
        // and not with BEAGLE's own rescaling. Otherwise an all-unknown tip stands in for the ghost subtree, which is
        // exact too, and keeps whatever resource the user chose.
        useOneChildOperations = IS_MATRIX_VECTOR_SUPPORTED()
                && (instanceFlags & BeagleFlag.FRAMEWORK_CPU.getMask()) != 0
                && !useAutoScaling
                && Boolean.parseBoolean(System.getProperty(GHOST_ONE_CHILD_PROPERTY, "true"));

        final String oneChildMethod = useOneChildOperations ? "from its other child alone" :
                "with an all-unknown tip in place of the ghost";
        if (mergeGhostBranches) {
            logger.info("    " + countGhosts(isGhostTip) + " of the " + tipCount + " taxa in the tree have no " +
                    "sequences in this partition; subtrees of them alone are skipped, and a node with one ghost child " +
                    "is skipped too, its branch joined to its other child's (the root, if it has a ghost child, " +
                    "is computed " + oneChildMethod + ").");
        } else {
            logger.info("    " + countGhosts(isGhostTip) + " of the " + tipCount + " taxa in the tree have no " +
                    "sequences in this partition; subtrees of them alone are skipped, and a node with one ghost child " +
                    "is computed " + oneChildMethod + ". Its branches cannot be joined because " + mergeProblem + ".");
        }
    }

    /**
     * @return why taxa of the tree that are missing from a partition (ghosts) cannot be used with these settings, as
     * the end of a sentence, or null if they can
     */
    public static String ghostProblem(PreOrderSettings settings) {
        if (settings.usePreOrder || settings.branchInfinitesimalDerivative) {
            // these read the partials, matrices and scale factors of every node, and ghosts have none
            return "that cannot be combined with pre-order partials or gradients (usePreOrder=\"true\")";
        }
        if (settings.useAugmentedEpochNodes) {
            return "that cannot yet be combined with augmented epoch nodes (useAugmentedEpochNodes=\"true\")";
        }
        return null;
    }

    /**
     * Joining the two branches through a node with one ghost child into one branch, of the summed length, is exact
     * only if every branch has the same rate matrix, so that a branch's transition matrix depends on its length
     * alone.
     *
     * @return why the branches through nodes with one ghost child cannot be joined, as the end of a sentence, or
     * null if they can
     */
    public static String ghostMergeProblem(BranchModel branchModel, PreOrderSettings settings) {
        if (!Boolean.parseBoolean(System.getProperty(GHOST_MERGE_PROPERTY, "true"))) {
            return "the property " + GHOST_MERGE_PROPERTY + " is false";
        }
        if (settings.useRewardAwareBranchModelDelegate) {
            return "the reward-aware branch model delegate is in use";
        }
        if (branchModel instanceof EpochBranchModel) {
            return "the branch model, " + branchModel.getId() + ", is an epoch model, with a different rate matrix " +
                    "in each epoch";
        }
        if (branchModel.getSubstitutionModels().size() != 1) {
            return "the branch model, " + branchModel.getId() + ", has more than one substitution model";
        }
        return null;
    }

    /**
     * Sets the status of each node that is being updated from its children's, which come before it.
     *
     * @return how many nodes changed from GHOST to having data, now listed at the start of becameReal
     */
    private int updateGhostStatus(List<NodeOperation> nodeOperations, int rootNodeNumber) {
        leavingScaledCount = 0;
        becameRealCount = 0;

        for (NodeOperation op : nodeOperations) {
            final int nodeNum = op.getNodeNumber();
            final boolean leftGhost = nodeStatus[op.getLeftChild()] == GHOST;
            final boolean rightGhost = nodeStatus[op.getRightChild()] == GHOST;
            final byte status;
            if (leftGhost == rightGhost) {
                status = leftGhost ? GHOST : FULL;
            } else {
                status = (mergeGhostBranches && nodeNum != rootNodeNumber) ? PASS : ONE_CHILD;
                if (mergeGhostBranches) {
                    passChild[nodeNum] = leftGhost ? op.getRightChild() : op.getLeftChild();
                }
            }
            setGhostStatus(nodeNum, status);
        }

        if (mergeGhostBranches) {
            // Whether a node with one ghost child is computed depends on whether it is the root, so a change of root
            // can change the status of a node whose children have not changed.
            if (previousRoot >= 0 && previousRoot != rootNodeNumber && nodeStatus[previousRoot] == ONE_CHILD) {
                setGhostStatus(previousRoot, PASS);
            }
            if (nodeStatus[rootNodeNumber] == PASS) {
                setGhostStatus(rootNodeNumber, ONE_CHILD);
            }
            previousRoot = rootNodeNumber;
        }

        return becameRealCount;
    }

    private void setGhostStatus(int nodeNum, byte status) {
        final byte oldStatus = nodeStatus[nodeNum];
        if (status == oldStatus) {
            return;
        }
        if (oldStatus == GHOST) {
            becameReal[becameRealCount++] = nodeNum;
        }

        // A node that stops being computed takes its factor out of the total at the root. One that stays computed
        // keeps it, as a node keeps its factor through any other change below it.
        if (hasScaleFactor[nodeNum] && !isComputed(status)) {
            leavingScaled[leavingScaledCount++] = scaleBufferIndices[nodeNum - tipCount];
            hasScaleFactor[nodeNum] = false;
        }

        if (mergeGhostBranches) {
            // its matrix buffer, if it had one in use, may be for another branch or out of date
            mergedLength[nodeNum] = Double.NaN;
            statusChangedAt[nodeNum] = evaluationStamp;
        }
        nodeStatus[nodeNum] = status;
    }

    /**
     * @return whether a node of this status is computed: every node but a ghost or a node passed through
     */
    private static boolean isComputed(byte status) {
        return status == FULL || status == ONE_CHILD;
    }

    /**
     * When the scale factors are not being recomputed, the total of them at the root is kept from the last time they
     * were. Scale factors are only written when they are recomputed, as without ghosts; between times a node keeps
     * its factor, wherever it moves in the tree, until it stops being computed, when its factor is taken out of the
     * total here. This happens on most moves once there are many ghosts, so the total is corrected rather than
     * recomputed. It is written to its other buffer, so that the stored one is still there to restore.
     */
    private void removeOldScaleFactors() {
        final int previous = scaleBufferHelper.getOffsetIndex(internalNodeCount);
        scaleBufferHelper.flipOffset(internalNodeCount);
        final int cumulative = scaleBufferHelper.getOffsetIndex(internalNodeCount);
        if (cumulative != previous) {
            beagle.copyScaleFactors(cumulative, previous);
        }
        if (leavingScaledCount > 0) {
            beagle.removeScaleFactors(leavingScaled, leavingScaledCount, cumulative);
        }
    }

    /**
     * Lists the branches whose matrices are needed: those of the updated nodes that have data, and those of nodes
     * that have just stopped being ghosts, whose matrices were not kept up to date while they were. The length of
     * every branch is remembered, ghost or not, for when it is needed.
     *
     * @return the number of branches listed in branchUpdateIndices and branchLengths
     */
    private int collectRealBranches(List<BranchOperation> branchOperations, int becameRealCount, int rootNodeNumber) {
        int count = 0;
        for (BranchOperation op : branchOperations) {
            final int branch = op.getBranchNumber();
            ghostBranchLengths[branch] = op.getBranchLength();
            if (nodeStatus[branch] != GHOST) {
                branchUpdateIndices[count] = branch;
                branchLengths[count] = op.getBranchLength();
                branchListed[branch] = true;
                count++;
            }
        }

        for (int i = 0; i < becameRealCount; i++) {
            final int node = becameReal[i];
            if (node != rootNodeNumber && !branchListed[node]) {
                branchUpdateIndices[count] = node;
                branchLengths[count] = ghostBranchLengths[node];
                count++;
            }
        }

        for (int i = 0; i < count; i++) {
            branchListed[branchUpdateIndices[i]] = false;
        }
        return count;
    }

    /**
     * Finds, for each node computed in this evaluation, its nearest descendant with data on each side, and lists the
     * joined branches whose matrices are needed. A node with data has one matrix buffer, its own, which holds the
     * matrix of the joined branch from it up to its nearest computed ancestor: its own branch and those of the nodes
     * passed through above it, which, under one rate matrix, is the matrix of their summed length.
     * <p>
     * A joined branch is listed if its length differs from the one its buffer was made for, or if the substitution
     * model or the category rates changed. Every change to a branch on it brings its top node into the node
     * operations, so nothing is missed.
     * <p>
     * It also decides which of these nodes are recomputed: those whose status, effective children, their partials or
     * the matrices of their joined branches changed. The others keep their partials.
     *
     * @return the number of branches listed in branchUpdateIndices and branchLengths
     */
    private int collectMergedBranches(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations) {
        for (BranchOperation op : branchOperations) {
            ghostBranchLengths[op.getBranchNumber()] = op.getBranchLength();
        }

        // A matrix depends only on the length of its branch, the substitution model and the rates of the site
        // categories, so a joined branch keeps its matrix unless one of these changed. (A branch operation alone is
        // not a change: one is sent for every child of a node whose children changed.) The category rates are
        // compared with what BEAGLE has, as every restore flags the site model as changed.
        if (updateSubstitutionModel) {
            allMatricesChanged = true;
        } else if (updateSiteModel) {
            final double[] categoryRates = siteRateModel.getCategoryRates();
            allMatricesChanged = categoryRates == null || !Arrays.equals(categoryRates, lastCategoryRates);
        } else {
            allMatricesChanged = false;
        }

        // every computed node is recomputed at a full rescale, and when the substitution model changed (as it has
        // after makeDirty, so that a full evaluation is one)
        final boolean recomputeAll = isFullRescale() || updateSubstitutionModel;

        int count = 0;
        for (NodeOperation op : nodeOperations) {
            final int nodeNum = op.getNodeNumber();
            final byte status = nodeStatus[nodeNum];
            if (status != FULL && status != ONE_CHILD) {
                continue;
            }
            final int oldLeft = effectiveLeft[nodeNum];
            final int oldRight = effectiveRight[nodeNum];
            if (status == FULL) {
                effectiveLeft[nodeNum] = op.getLeftChild();
                effectiveRight[nodeNum] = op.getRightChild();
            } else {
                effectiveLeft[nodeNum] = passChild[nodeNum];
                effectiveRight[nodeNum] = NodeOperation.NO_CHILD;
            }
            final int before = count;
            count = listJoinedBranch(effectiveLeft, nodeNum, count);
            if (status == FULL) {
                count = listJoinedBranch(effectiveRight, nodeNum, count);
            }
            final int left = effectiveLeft[nodeNum];
            final int right = effectiveRight[nodeNum];

            // A node's partials depend only on its status, its effective children, their partials and the matrices
            // of their joined branches. If none of these changed, which is so when a move touched ghosts alone, its
            // partials would come out the same, so it keeps them.
            final boolean sameChildren = (left == oldLeft && right == oldRight) || (left == oldRight && right == oldLeft);
            if (recomputeAll || statusChangedAt[nodeNum] == evaluationStamp || !sameChildren || count > before ||
                    isRecomputed(left) || isRecomputed(right)) {
                recomputedAt[nodeNum] = evaluationStamp;
            }
        }

        return count;
    }

    /**
     * Walks down from the child in effective[top] through any nodes passed through to the first node with data on
     * both sides, or a tip, leaving that node in effective[top], and lists its joined branch if it needs a new matrix.
     */
    /**
     * @return whether every computed node must write a new scale factor: a periodic rescale or the retry after an
     * underflow. With rescaling ALWAYS, factors are recomputed on every evaluation, but only for the nodes that are
     * computed: a node whose partials do not change keeps the factor that goes with them.
     */
    private boolean isFullRescale() {
        return recomputeScaleFactors && rescalingScheme == PartialsRescalingScheme.DYNAMIC;
    }

    /**
     * @return whether a node's partials are being recomputed in this evaluation (never so for a tip or NO_CHILD)
     */
    private boolean isRecomputed(int node) {
        return node >= 0 && recomputedAt[node] == evaluationStamp;
    }

    private int listJoinedBranch(int[] effective, int top, int count) {
        int node = effective[top];
        double length = ghostBranchLengths[node];
        while (nodeStatus[node] == PASS) {
            node = passChild[node];
            length += ghostBranchLengths[node];
        }
        effective[top] = node;

        if (allMatricesChanged || length != mergedLength[node]) { // a NaN length is never equal, so is always listed
            mergedLength[node] = length;
            branchUpdateIndices[count] = node;
            branchLengths[count] = length;
            count++;
        }
        return count;
    }

    /**
     * Lists the scale buffers of the nodes whose partials carry a scale factor of their own: with BEAGLE's own
     * scaling, every node that is computed; otherwise those that have had one since the factors were last recomputed.
     *
     * @return the number of buffers listed in scaledNodeBufferIndices
     */
    private int collectScaledNodes() {
        int count = 0;
        for (int i = 0; i < internalNodeCount; i++) {
            if (useAutoScaling ? isComputed(nodeStatus[tipCount + i]) : hasScaleFactor[tipCount + i]) {
                scaledNodeBufferIndices[count++] = scaleBufferIndices[i];
            }
        }
        return count;
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    // The status of a node in a partition with ghosts: FULL has data on both sides, ONE_CHILD on one side only, and
    // a GHOST node has none below it at all. FULL is zero so that new arrays start with it.
    private static final byte FULL = 0;
    private static final byte ONE_CHILD = 1;
    private static final byte GHOST = 2;
    // A node with one ghost child, other than the root, when the branches through it are joined: it is not computed,
    // and its parent reads its real child's partials across the joined branch. A ONE_CHILD node is then the root.
    private static final byte PASS = 3;

    // null when every taxon of the tree is in the pattern list, which then takes the ordinary path throughout
    private final boolean[] isGhostTip;
    private byte[] nodeStatus;
    private byte[] storedNodeStatus;
    private double[] ghostBranchLengths;        // the length of every branch, for nodes that stop being ghosts
                                                // and for summing the branches through nodes passed through
    private double[] storedGhostBranchLengths;
    private int[] becameReal;                   // scratch: nodes that stopped being ghosts in this evaluation
    private int becameRealCount;
    private boolean[] branchListed;             // scratch: branches already listed for a matrix update
    private int[] scaledNodeBufferIndices;      // scratch: the scale buffers to accumulate at the root
    private boolean[] hasScaleFactor;           // whether a node's partials are divided by a factor in the total
    private boolean[] storedHasScaleFactor;
    private int[] leavingScaled;                // scratch: the scale buffers of nodes that stopped being computed
    private int leavingScaledCount;
    private boolean useOneChildOperations;
    private int ghostTipNumber = -1;            // a ghost tip, all unknown, that can stand in for a ghost subtree

    // Joining the branches through nodes with one ghost child. The matrix buffer of a node with data then holds the
    // matrix of the joined branch from it up to its nearest ancestor that is computed (FULL, or the root).
    private boolean mergeGhostBranches;
    private int[] passChild;                    // for a PASS or ONE_CHILD node, its child with data
    private int[] storedPassChild;
    private double[] mergedLength;              // the length a node's matrix buffer was made for, or NaN
    private double[] storedMergedLength;
    private int[] effectiveLeft;                // the nearest descendants with data of a computed node
    private int[] effectiveRight;
    private boolean allMatricesChanged;         // scratch: whether the model that every matrix depends on changed
    private int previousRoot = -1;
    private int storedPreviousRoot = -1;
    private int[] storedEffectiveLeft;
    private int[] storedEffectiveRight;
    private int evaluationStamp = 0;            // counts the evaluations, to mark nodes without clearing arrays
    private int[] recomputedAt;                 // the evaluation in which a computed node was last recomputed
    private int[] statusChangedAt;              // the evaluation in which a node's status last changed
    private double lastLogL = Double.NaN;       // the log likelihood of the last evaluation that computed one
    private double storedLastLogL = Double.NaN;
    private int lastRootNode = -1;              // and the root it was computed at
    private int storedLastRootNode = -1;
    private boolean rootInputsChanged;          // whether the category rates or weights or root frequencies changed
    private double[] lastCategoryRates;         // what the matrices and the last likelihood were computed with
    private double[] lastCategoryWeights;
    private double[] lastRootFrequencies;
    private double[] storedLastCategoryRates;
    private double[] storedLastCategoryWeights;
    private double[] storedLastRootFrequencies;

    private final int nodeCount;
    private final int tipCount;
    private final int internalNodeCount;

    private final int[] branchUpdateIndices;
    private final double[] branchLengths;

    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private final int[] operations;

    private EpochTimeProvider epochTimeProvider = null; // non-null when there are augmented epoch nodes

    // non-null when there are augmented epoch nodes
    private final EpochSubstitutionModelDelegate epochDelegate;
    private final AugmentedNodeRegistry registry;
    private final AugmentedBufferIndexHelper augmentedPartialHelper; // the same object as partialBufferHelper

    private final boolean growable;  // the BEAGLE buffers of augmented nodes grow on demand (BEAGLE 4.2.0)
    private int allocatedAugmented;  // augmented nodes with BEAGLE buffers
    private int growthEventCount;

    private final int compactPartialsCount;

    private boolean flip = true;
    private final BufferIndexHelper partialBufferHelper;
    private final BufferIndexHelper scaleBufferHelper;

    private PartialsRescalingScheme rescalingScheme;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private boolean delayRescalingUntilUnderflow = true;

    private boolean useScaleFactors = false;
    private boolean useAutoScaling = false;

    private boolean recomputeScaleFactors = false;
    // incremented by each likelihood evaluation that writes individual scale factors
    private long scaleFactorsVersion = 0;
    // per gradient scale buffer slot (isNegAncestorScaleBufferCurrent): its source and individual buffers, the version
    // and the source's build it was made from, and its own build (a count over all slots; 0: never built)
    private int[] negAncestorSource = new int[0];
    private int[] negAncestorIndividual = new int[0];
    private long[] negAncestorVersion = new long[0];
    private long[] negAncestorSourceBuild = new long[0];
    private long[] negAncestorBuild = new long[0];
    private long negAncestorBuildCount = 0;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    private int threadCount = -1;
    private long instanceFlags;

    private boolean firstRescaleAttempt = false;
    private int rescalingMessageCount = 0;

    //integer to keep track of setting recomputeScaleFactors correctly after an underflow
    private int underflowHandling = 0;

    /**
     * the patternList
     */
    private final PatternList patternList;

    /**
     * the data type
     */
    private final DataType dataType;

    /**
     * the pattern weights
     */
    private double[] patternWeights;
    private double[] storedPatternWeights;
    private Parameter siteAssignInd;
    private Parameter partitionCat;

    /**
     * if the chain is already restored
     */
    private boolean isRestored;

    /**
     * the number of patterns
     */
    private final int patternCount;

    /**
     * the number of states in the data
     */
    private final int stateCount;

    /**
     * the branch-site model for these sites
     */
    private final BranchModel branchModel;

    /**
     * A delegate to handle substitution models on branches
     */
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;

    /**
     * the site model for these sites
     */
    private final SiteRateModel siteRateModel;

    /**
     * the pattern likelihoods
     */
    private double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    private final int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    private double[] tipPartials;

    /**
     * an array used to transfer tip states
     */
    private int[] tipStates;

    /**
     * the BEAGLE library instance
     */
    private final Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    private boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    private boolean updateSiteModel;

    /**
     * Flag to specify that the root frequencies has changed
     */
    private boolean updateRootFrequency;

    /**
     * Flag to take into account the first likelihood evaluation when initiating the MCMC chain
     */
    private boolean initialEvaluation = true;

    /**
     * PreOrder related settings
     */
    private PreOrderSettings settings;

    private boolean useAmbiguities;
    private boolean prefersGPU;
    private static boolean USE_CACHED_EXCEPTION = true;
    private LikelihoodUnderflowException cachedException = null; // new LikelihoodUnderflowException();

    @Override
    public void setTipStates(int tipNum, int[] states) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void getTipStates(int tipNum, int[] states) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getPatternCount() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getTipCount() {
        return tipCount;
    }

    // **************************************************************
    // Threading (used to share a chain's processor budget between its instances)
    // **************************************************************

    /**
     * @return the number of unique site patterns computed by this instance
     */
    public int getInstancePatternCount() {
        return patternCount;
    }

    /**
     * @return the number of states of the data computed by this instance
     */
    public int getInstanceStateCount() {
        return stateCount;
    }

    /**
     * @return true if this instance runs on a CPU (rather than a GPU or other accelerator)
     */
    public boolean isCPUInstance() {
        return (instanceFlags & BeagleFlag.PROCESSOR_CPU.getMask()) != 0;
    }

    /**
     * @return the number of CPU threads this instance has been asked to use (-1 for BEAGLE's
     * own choice, 0 or 1 for none)
     */
    public int getCPUThreadCount() {
        return threadCount;
    }

    /**
     * Changes the number of threads used by a CPU instance. This only has an effect if the
     * instance was created with threading on.
     *
     * @return true if the thread count was changed
     */
    public boolean setCPUThreadCount(int threadCount) {
        if (isCPUInstance() && IS_THREAD_COUNT_COMPATIBLE() &&
                (instanceFlags & BeagleFlag.THREADING_CPP.getMask()) != 0 && threadCount > 0) {
            beagle.setCPUThreadCount(threadCount);
            this.threadCount = threadCount;
            return true;
        }
        return false;
    }
}
