/*
 * BeagleDataLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood;

import beagle.*;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodel.treedatalikelihood.BeagleFunctionality.*;

/**
 * BeagleDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate that uses BEAGLE
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

public class BeagleDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate, Citable {

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

    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

    public static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;
    private static List<Integer> extraBufferOrder = null;

    // Default frequency for complete recomputation of scaling factors under the 'dynamic' scheme
    private static final int RESCALE_FREQUENCY = 100;
    private static final int RESCALE_TIMES = 1;

    // count the number of partial likelihood and matrix updates
    private long totalMatrixUpdateCount = 0;
    private long totalPartialsUpdateCount = 0;
    private long totalEvaluationCount = 0;

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
                                        PreOrderSettings settings) {

        super("BeagleDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing BEAGLE DataLikelihood Delegate");
        setId(patternList.getId());

        this.dataType = patternList.getDataType();
        this.patternList = patternList;
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

        // Check for matching state counts
        int stateCount2 = branchModel.getRootFrequencyModel().getFrequencyCount();
        if (stateCount != stateCount2) {
            throw new IllegalArgumentException("Pattern state count (" + stateCount
                    + ") does not match substitution model state count (" + stateCount2 + ")");
        }

        patternWeights = patternList.getPatternWeights();

        this.branchModel = branchModel;
        addModel(this.branchModel);

        this.siteRateModel = siteRateModel;
        addModel(this.siteRateModel);

        this.categoryCount = this.siteRateModel.getCategoryCount();

        nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        branchUpdateIndices = new int[nodeCount];
        branchLengths = new double[nodeCount];
        scaleBufferIndices = new int[internalNodeCount];
        storedScaleBufferIndices = new int[internalNodeCount];

        operations = new int[internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];

        firstRescaleAttempt = true;

        isRestored = false;

        this.settings = settings;

        try {

            int compactPartialsCount = tipCount;
            if (useAmbiguities) {
                // if we are using ambiguities then we don't use tip partials
                compactPartialsCount = 0;
            }

            // one partials buffer for each tip and two for each internal node (for store restore)
            partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

            // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
            scaleBufferHelper = new BufferIndexHelper(getSingleScaleBufferCount(), 0);

            if (settings.branchInfinitesimalDerivative) {
                evolutionaryProcessDelegate = new SubstitutionModelDelegate(tree, branchModel, settings);
            } else {

                if (branchModel.getSubstitutionModels().size() == 1) {
                    evolutionaryProcessDelegate = new HomogenousSubstitutionModelDelegate(tree, branchModel);
                } else {
                    // use a more general delegate that allows different substitution models on different branches and
                    // can do matrix convolution.

                    // TODO: the constructor should take the delegate and the delegate should wrap the branchModel
                    evolutionaryProcessDelegate = new SubstitutionModelDelegate(tree, branchModel, settings);
                }
            }

            int numPartials = partialBufferHelper.getBufferCount();
            int numScaleBuffers = scaleBufferHelper.getBufferCount();
            int numMatrices = evolutionaryProcessDelegate.getMatrixBufferCount();

            // one partial buffer for root node and two for each node including tip nodes (for store restore)
            if (settings.usePreOrder){
                numPartials += nodeCount;
                numScaleBuffers += nodeCount - 1; // don't need to rescale at root
                numMatrices += evolutionaryProcessDelegate.getCachedMatrixBufferCount(settings);
            }

            // Attempt to get the resource order from the System Property
            if (resourceOrder == null) {
                resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
            }
            if (preferredOrder == null) {
                preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
            }
            if (requiredOrder == null) {
                requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
            }
            if (scalingOrder == null) {
                scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
            }
            if (extraBufferOrder == null) {
                extraBufferOrder = parseSystemPropertyIntegerArray(EXTRA_BUFFER_COUNT_PROPERTY);
            }

            // first set the rescaling scheme to use from the parser
            this.rescalingScheme = rescalingScheme;
            this.delayRescalingUntilUnderflow = delayRescalingUntilUnderflow;

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
                if (resourceList[0] > 0) {
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

            if (threadCount == 0 || threadCount == 1) {
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
                    (resourceList != null && resourceList[0] > 0)) {
                // non-CPU implementations don't have SSE so remove default preference for SSE
                // when using non-CPU preferences or prioritising non-CPU resource
                preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
                preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
            }

            if (settings.usePreOrder && stateCount > 4) {
                requirementFlags |= BeagleFlag.PREORDER_TRANSPOSE_AUTO.getMask();
            }

            // start auto resource selection
            String resourceAuto = System.getProperty(RESOURCE_AUTO_PROPERTY);
            if (resourceAuto != null && Boolean.parseBoolean(resourceAuto)) {

                long benchmarkFlags = 0;

                if (this.rescalingScheme == PartialsRescalingScheme.NONE) {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_NONE.getMask();
                } else if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS) {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_ALWAYS.getMask();
                } else {
                    benchmarkFlags =  BeagleBenchmarkFlag.SCALING_DYNAMIC.getMask();
                }

                logger.info("\nRunning benchmarks to automatically select fastest BEAGLE resource for analysis or partition... ");

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

            if (IS_THREAD_COUNT_COMPATIBLE() && threadCount > 1) {
                beagle.setCPUThreadCount(threadCount);
            }

            if (patternList instanceof UncertainSiteList) { // TODO Remove
                useAmbiguities = true;
            }

//            if (!IS_PRE_ORDER_SUPPORTED() && settings.usePreOrder) {
//                throw new IllegalArgumentException("BEAGLE library does not support pre-order computation");
//            }

            //add in logger info for preOrder traversal
            logger.info("  " + (settings.usePreOrder ? "Using" : "Ignoring") + " preOrder partials in tree likelihood.");
            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");

            if (patternList.areUncertain() && !useAmbiguities) {
                logger.info("  WARNING: Uncertain site patterns will be ignored.");
            }

            for (int i = 0; i < tipCount; i++) {
                // Find the id of tip i in the patternList
                String id = tree.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + tree.getId() +
                            ", is not found in patternList, " + patternList.getId());
                } else {
                    if (useAmbiguities) {
                        setPartials(beagle, patternList, index, i);
                    } else {
                        setStates(beagle, patternList, index, i);
                    }
                }
            }

            beagle.setPatternWeights(patternWeights);

            String rescaleMessage = "  Using rescaling scheme : " + this.rescalingScheme.getText();
            if (this.rescalingScheme == PartialsRescalingScheme.AUTO &&
                    resourceDetails != null &&
                    (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
                // If auto scaling in BEAGLE is not supported then do it here
                this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                rescaleMessage = "  Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText();
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

        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchLength();
            branchUpdateCount ++;
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
            // these could be set only when they change but store/restore would need to be considered
            beagle.setCategoryWeights(0, categoryWeights); // TODO move
        }

        if (updateRootFrequency) {
            double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
            beagle.setStateFrequencies(0, frequencies); // TODO make lazy?
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
                partialBufferHelper.flipOffset(op.getNodeNumber());
            }
        }

        int operationCount = nodeOperations.size();
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            int nodeNum = op.getNodeNumber();

            operations[k] = partialBufferHelper.getOffsetIndex(nodeNum);

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

            operations[k + 3] = partialBufferHelper.getOffsetIndex(op.getLeftChild()); // source node 1
            operations[k + 4] = evolutionaryProcessDelegate.getMatrixIndex(op.getLeftChild()); // source matrix 1
            operations[k + 5] = partialBufferHelper.getOffsetIndex(op.getRightChild()); // source node 2
            operations[k + 6] = evolutionaryProcessDelegate.getMatrixIndex(op.getRightChild()); // source matrix 2

            k += Beagle.OPERATION_TUPLE_SIZE;
        }

        beagle.updatePartials(operations, operationCount, Beagle.NONE);

        if (COUNT_CALCULATIONS) {
            totalEvaluationCount += 1;
            totalPartialsUpdateCount += operationCount;
        }

        int rootIndex = partialBufferHelper.getOffsetIndex(rootNodeNumber);

        int cumulateScaleBufferIndex = Beagle.NONE;

        if (useScaleFactors) {
            if (recomputeScaleFactors) {
                scaleBufferHelper.flipOffset(internalNodeCount);
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                beagle.resetScaleFactors(cumulateScaleBufferIndex);
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
            } else {
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
            }
        } else if (useAutoScaling) {
            beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, Beagle.NONE);
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
                            Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling... (" + getId() + ")");
                        }
                    }
                    rescalingMessageCount += 1;
                }

                useScaleFactors = true;
                recomputeScaleFactors = true;

                firstRescaleAttempt = false; // Only try to rescale once

                rescalingCount--;

            }

            // turn off double buffer flipping so the next call overwrites the
            // underflowed buffers. Flip will be turned on again in storeState for
            // next step
            flip = false;
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

    public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        /* No need to rescale partials */
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
    }

    private void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }

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
        } else if (model == branchModel) {
            updateSubstitutionModel = true;
            updateRootFrequency = true;
        }

        // Tell TreeDataLikelihood to update all nodes
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void storeState() {
        partialBufferHelper.storeState();
        evolutionaryProcessDelegate.storeState();

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
        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore
        updateRootFrequency = true;

        partialBufferHelper.restoreState();
        evolutionaryProcessDelegate.restoreState();

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

    public final EvolutionaryProcessDelegate getEvolutionaryProcessDelegate(){
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

    public final int getPartialBufferCount() {
        return partialBufferHelper.getBufferCount();
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
        if (beagle != null) {
            beagle.finalize();
        }
    }

    public static void releaseAllBeagleInstances() throws Throwable {
        for (Likelihood likelihood : dr.inference.model.Likelihood.FULL_LIKELIHOOD_SET) {
            if (likelihood instanceof TreeDataLikelihood) {
                TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) likelihood;
                DataLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getDataLikelihoodDelegate();
                if (likelihoodDelegate instanceof BeagleDataLikelihoodDelegate) {
                    BeagleDataLikelihoodDelegate delegate = (BeagleDataLikelihoodDelegate) likelihoodDelegate;
                    delegate.releaseBeagle();
                }
            }
        }
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final int nodeCount;
    private final int tipCount;
    private final int internalNodeCount;

    private final int[] branchUpdateIndices;
    private final double[] branchLengths;

    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private final int[] operations;

    private boolean flip = true;
    private final BufferIndexHelper partialBufferHelper;
    private final BufferIndexHelper scaleBufferHelper;

    private PartialsRescalingScheme rescalingScheme;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private boolean delayRescalingUntilUnderflow = true;

    private boolean useScaleFactors = false;
    private boolean useAutoScaling = false;

    private boolean recomputeScaleFactors = false;
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
    private final double[] patternWeights;

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

    private static boolean USE_CACHED_EXCEPTION = true;
    private LikelihoodUnderflowException cachedException = null; // new LikelihoodUnderflowException();
}
