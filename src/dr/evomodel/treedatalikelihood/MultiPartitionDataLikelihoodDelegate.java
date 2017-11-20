/*
 * MultiPartitionDataLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * MultiPartitionDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate that uses BEAGLE 3 to allow for parallelization across multiple data partitions
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @author Guy Baele
 * @version $Id$
 */

import beagle.*;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

public class MultiPartitionDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate, Citable {

    private static final boolean RESCALING_OFF = false; // a debugging switch
    private static final boolean DEBUG = false;

    public static boolean IS_MULTI_PARTITION_COMPATIBLE() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return versionNumbers.length != 0 && versionNumbers[0] >= 3;
    }

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    private static final String DELAY_SCALING_PROPERTY = "beagle.delay.scaling";
    private static final String EXTRA_BUFFER_COUNT_PROPERTY = "beagle.extra.buffer.count";
    private static final String FORCE_VECTORIZATION = "beagle.force.vectorization";

    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;
    private static List<Integer> extraBufferOrder = null;

    // Default frequency for complete recomputation of scaling factors under the 'dynamic' scheme
    private static final int RESCALE_FREQUENCY = 100;
    private static final int RESCALE_TIMES = 1;

    /**
     * Construct an instance using a list of PatternLists, one for each partition. The
     * partitions will share a tree but can have different branchModels and siteRateModels
     * The latter should either have a size of 1 (in which case they are shared across partitions)
     * or equal to patternLists.size() where each partition has a different model.
     *
     * @param tree Used for configuration - shouldn't be watched for changes
     * @param branchModels Specifies a list of branch models for each partition
     * @param patternLists List of patternLists comprising each partition
     * @param siteRateModels A list of siteRateModels for each partition
     * @param useAmbiguities Whether to respect state ambiguities in data
     */
    public MultiPartitionDataLikelihoodDelegate(Tree tree,
                                                List<PatternList> patternLists,
                                                List<BranchModel> branchModels,
                                                List<SiteRateModel> siteRateModels,
                                                boolean useAmbiguities,
                                                PartialsRescalingScheme rescalingScheme,
                                                boolean delayRescalingUntilUnderflow) {

        super("MultiPartitionDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing Multi-Partition Data Likelihood Delegate with BEAGLE 3 extensions");

        setId(patternLists.get(0).getId());

        this.patternLists = patternLists;
        this.dataType = patternLists.get(0).getDataType();
        stateCount = dataType.getStateCount();

        partitionCount = patternLists.size();
        patternCounts = new int[partitionCount];
        int total = 0;
        int k = 0;
        for (PatternList patternList : patternLists) {
            assert(patternList.getDataType().equals(this.dataType));
            patternCounts[k] = patternList.getPatternCount();
            total += patternCounts[k];
            k++;
        }
        totalPatternCount = total;

        useScaleFactors = new boolean[partitionCount];
        recomputeScaleFactors = new boolean[partitionCount];
        everUnderflowed = new boolean[partitionCount];
        flip = new boolean[partitionCount];
        for (int i = 0; i < partitionCount; i++) {
            flip[i] = true;
        }

        updatePartition = new boolean[partitionCount];
        partitionWasUpdated = new boolean[partitionCount];
        updateAllPartitions = true;

        cachedLogLikelihoodsByPartition = new double[partitionCount];
        storedCachedLogLikelihoodsByPartition = new double[partitionCount];

        // Branch models determine the substitution models per branch. There can be either
        // one per partition or one shared across all partitions
        assert(branchModels.size() == 1 || (branchModels.size() == patternLists.size()));

        this.branchModels.addAll(branchModels);
        for (BranchModel branchModel : this.branchModels) {
            addModel(branchModel);
        }

        // SiteRateModels determine the rates per category (for site-heterogeneity models).
        // There can be either one per partition or one shared across all partitions
        assert(siteRateModels.size() == 1 || (siteRateModels.size() == patternLists.size()));

        this.siteRateModels.addAll(siteRateModels);
        this.categoryCount = this.siteRateModels.get(0).getCategoryCount();
        for (SiteRateModel siteRateModel : this.siteRateModels) {
            assert(siteRateModel.getCategoryCount() == categoryCount);
            addModel(siteRateModel);
        }

        nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        branchUpdateIndices = new int[nodeCount];
        branchLengths = new double[nodeCount];
        //changed initialization to account for multiple partitions
        scaleBufferIndices = new int[partitionCount][internalNodeCount];
        storedScaleBufferIndices = new int[partitionCount][internalNodeCount];

        operations = new int[internalNodeCount * Beagle.PARTITION_OPERATION_TUPLE_SIZE * partitionCount];

        rescalingCount = new int[partitionCount];
        rescalingCountInner = new int[partitionCount];
        firstRescaleAttempt = true;

        try {

            int compactPartialsCount = tipCount;
            if (useAmbiguities) {
                // if we are using ambiguities then we don't use tip partials
                compactPartialsCount = 0;
            }

            partialBufferHelper = new BufferIndexHelper[partitionCount];
            scaleBufferHelper = new BufferIndexHelper[partitionCount];
            categoryRateBufferHelper = new BufferIndexHelper[partitionCount];
            for (int i = 0; i < partitionCount; i++) {
                // one partials buffer for each tip and two for each internal node (for store restore)
                partialBufferHelper[i] = new BufferIndexHelper(nodeCount, tipCount);

                // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
                scaleBufferHelper[i] = new BufferIndexHelper(getScaleBufferCount(), 0);

                categoryRateBufferHelper[i] = new BufferIndexHelper(1, 0);
            }

            int eigenBufferCount = 0;
            int matrixBufferCount = 0;

            // create a substitutionModelDelegate for each branchModel
            int partitionNumber = 0;
            for (BranchModel branchModel : this.branchModels) {
                HomogenousSubstitutionModelDelegate substitutionModelDelegate = new HomogenousSubstitutionModelDelegate(tree, branchModel, partitionNumber);
                evolutionaryProcessDelegates.add(substitutionModelDelegate);

                eigenBufferCount += substitutionModelDelegate.getEigenBufferCount();
                matrixBufferCount += substitutionModelDelegate.getMatrixBufferCount();

                partitionNumber ++;
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
                } else {
                    // if CPU: just run as fast as possible
                    // this.rescalingScheme = PartialsRescalingScheme.NONE;
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
                // auto scaling not supported for multi-partition BEAGLE3 instances
                preferenceFlags |= BeagleFlag.SCALING_DYNAMIC.getMask();
            } else {
                // preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
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

            // I don't think this performance stuff should be here. Perhaps have an intelligent automatic
            // load balancer further up the chain.
//            if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
//                if (stateCount == 4 && patternList.getPatternCount() < 10000) // TODO determine good cut-off
//                    preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
//            }

            boolean forceVectorization = false;
            String vectorizationString = System.getProperty(FORCE_VECTORIZATION);
            if (vectorizationString != null) {
                forceVectorization = true;
            }

            if (BeagleFlag.VECTOR_SSE.isSet(preferenceFlags) && (stateCount != 4)
                    && !forceVectorization
                    ) {
                // @todo SSE doesn't seem to work for larger state spaces so for now we override the
                // SSE option.
                preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
                preferenceFlags |= BeagleFlag.VECTOR_NONE.getMask();

                if (stateCount > 4 && this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                    this.rescalingScheme = PartialsRescalingScheme.DELAYED;
                }
            }

            if (!BeagleFlag.PRECISION_SINGLE.isSet(preferenceFlags)) {
                // if single precision not explicitly set then prefer double
                preferenceFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
            }

            if (evolutionaryProcessDelegates.get(0).canReturnComplexDiagonalization()) {
                requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
            }

            //TODO: check getBufferCount() calls with Daniel
            //TODO: should we multiple getBufferCount() by the number of partitions?
            beagle = BeagleFactory.loadBeagleInstance(
                    tipCount,
                    partialBufferHelper[0].getBufferCount(),
                    compactPartialsCount,
                    stateCount,
                    totalPatternCount,
                    eigenBufferCount,
                    matrixBufferCount,
                    categoryCount,
                    scaleBufferHelper[0].getBufferCount(), // Always allocate; they may become necessary
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

            patternPartitions = new int[totalPatternCount];
            patternWeights = new double[totalPatternCount];

            int j = 0;
            k = 0;
            for (PatternList patternList : patternLists) {
                double[] pw = patternList.getPatternWeights();
                for (int i = 0; i < patternList.getPatternCount(); i++) {
                    patternPartitions[k] = j;
                    patternWeights[k] = pw[i];
                    k++;
                }
                j++;
            }

            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            String patternCountString = "" + patternLists.get(0).getPatternCount();
            for (int i = 1; i < patternLists.size(); i++) {
                patternCountString += ", " + patternLists.get(i).getPatternCount();
            }
            logger.info("  With " + patternLists.size() + " partitions comprising " + patternCountString + " unique site patterns");

            // @todo - should check that each patternList spans the same set of taxa
            for (int i = 0; i < tipCount; i++) {
                String id = tree.getTaxonId(i);
                if (useAmbiguities) {
                    setPartials(beagle, patternLists, id, i);
                } else {
                    setStates(beagle, patternLists, id, i);
                }
            }

            beagle.setPatternWeights(patternWeights);
            beagle.setPatternPartitions(partitionCount, patternPartitions);

            String rescaleMessage = "  Using rescaling scheme : " + this.rescalingScheme.getText();
            if (this.rescalingScheme == PartialsRescalingScheme.AUTO) {
                // auto scaling in BEAGLE3 is not supported
                this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                rescaleMessage = "  Auto rescaling not supported in BEAGLE v3, using : " + this.rescalingScheme.getText();
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
                for (int i = 0; i < partitionCount; i++) {
                    everUnderflowed[i] = false; // If false, BEAST does not rescale until first under-/over-flow.
                }
            }

            updateSubstitutionModels = new boolean[branchModels.size()];
            updateSubstitutionModels();

            updateSiteRateModels = new boolean[siteRateModels.size()];
            updateSiteRateModels();

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
    }

    @Override
    public String getReport() {
        return null;
    }

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        return TreeTraversal.TraversalType.REVERSE_LEVEL_ORDER;
    }

    public List<PatternList> getPatternLists() {
        return this.patternLists;
    }

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getTraitDim() {
        return totalPatternCount;
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        return RateRescalingScheme.NONE;
    }

    private void updateSubstitutionModels(boolean... state) {
        for (int i = 0; i < updateSubstitutionModels.length; i++) {
            updateSubstitutionModels[i] = (state.length < 1 || state[0]);
        }
    }

    private void updateSubstitutionModel(BranchModel branchModel) {
        for (int i = 0; i < branchModels.size(); i++) {
            if (branchModels.get(i) == branchModel) {
                updateSubstitutionModels[i] = true;
            }
        }
    }

    private void updateSiteRateModels(boolean... state) {
        for (int i = 0; i < updateSiteRateModels.length; i++) {
            updateSiteRateModels[i] = (state.length < 1 || state[0]);
        }
    }

    private void updateSiteRateModel(SiteRateModel siteRateModel) {
        for (int i = 0; i < siteRateModels.size(); i++) {
            if (siteRateModels.get(i) == siteRateModel) {
                updateSiteRateModels[i] = true;
            }
        }
    }


    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<Integer>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<String>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private int getScaleBufferCount() {
        return internalNodeCount + 1;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternLists  patternLists
     * @param taxonId       taxonId
     * @param nodeIndex     nodeIndex
     */
    private final void setPartials(Beagle beagle,
                                   List<PatternList> patternLists,
                                   String taxonId,
                                   int nodeIndex) throws TaxonList.MissingTaxonException {

        double[] partials = new double[totalPatternCount * stateCount * categoryCount];
        int v = 0;
        for (PatternList patternList : patternLists) {
            int sequenceIndex = patternList.getTaxonIndex(taxonId);

            if (sequenceIndex == -1) {
                throw new TaxonList.MissingTaxonException("Taxon, " + taxonId +
                        ", not found in patternList, " + patternList.getId());
            }

            boolean[] stateSet;

            for (int i = 0; i < patternList.getPatternCount(); i++) {

                int state = patternList.getPatternState(sequenceIndex, i);
                stateSet = dataType.getStateSet(state);

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
        int n = totalPatternCount * stateCount;
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
     * @param patternLists  patternLists
     * @param taxonId       taxonId
     * @param nodeIndex     nodeIndex
     */
    private final void setStates(Beagle beagle,
                                 List<PatternList> patternLists,
                                 String taxonId,
                                 int nodeIndex) throws TaxonList.MissingTaxonException {

        int[] states = new int[totalPatternCount];

        int v = 0;
        for (PatternList patternList : patternLists) {
            int sequenceIndex = patternList.getTaxonIndex(taxonId);

            if (sequenceIndex == -1) {
                throw new TaxonList.MissingTaxonException("Taxon, " + taxonId +
                        ", not found in patternList, " + patternList.getId());
            }

            for (int i = 0; i < patternList.getPatternCount(); i++) {
                states[v] = patternList.getPatternState(sequenceIndex, i);
                v++;
            }
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

        boolean throwLikelihoodRescalingException = false;
        if (!initialEvaluation) {
            for (int i = 0; i < partitionCount; i++) {
                if (!this.delayRescalingUntilUnderflow || everUnderflowed[i]) {
                    if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS || this.rescalingScheme == PartialsRescalingScheme.DELAYED) {
                        useScaleFactors[i] = true;
                        recomputeScaleFactors[i] = true;
                    } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                        useScaleFactors[i] = true;

                        if (DEBUG) {
                            System.out.println("rescalingCount["+i+"] = " + rescalingCount[i]);
                        }

                        if (rescalingCount[i] > rescalingFrequency) {
                            if (DEBUG) {
                                System.out.println("rescalingCount > rescalingFrequency");
                            }
                            rescalingCount[i] = 0;
                            rescalingCountInner[i] = 0;
                        }

                        if (DEBUG) {
                            System.out.println("rescalingCountInner = " + rescalingCountInner[i]);
                        }

                        if (rescalingCountInner[i] < RESCALE_TIMES) {
                            if (DEBUG) {
                                System.out.println("rescalingCountInner < RESCALE_TIMES");
                            }

                            recomputeScaleFactors[i] = true;
                            updatePartition[i] = true;

                            rescalingCountInner[i]++;

                            throwLikelihoodRescalingException = true;
                        }
                    }
                }
            }

            if (throwLikelihoodRescalingException) {
                throw new LikelihoodRescalingException();
            }
        }

        if (RESCALING_OFF) { // a debugging switch
            for (int i = 0; i < partitionCount; i++) {
                useScaleFactors[i] = false;
                recomputeScaleFactors[i] = false;
            }
        }

        int k = 0;
        for (EvolutionaryProcessDelegate evolutionaryProcessDelegate : evolutionaryProcessDelegates) {
            if (updateSubstitutionModels[k]) {
                // TODO: More efficient to update only the substitution model that changed, instead of all
                // TODO: flip currently assumes 1 substitution model per partition
                evolutionaryProcessDelegate.updateSubstitutionModels(beagle, flip[k]);
                updatePartition[k] = true;
                if (DEBUG) {
                    System.out.println("updateSubstitutionModels, updatePartition["+k+"] = " + updatePartition[k]);
                }
                updateAllPartitions = false;
                // we are currently assuming a no-category model...
            }
            k++;
        }

        k = 0;
        for (SiteRateModel siteRateModel : siteRateModels) {
            if (updateSiteRateModels[k]) {
                double[] categoryRates = siteRateModel.getCategoryRates();
                if (categoryRates == null) {
                    // If this returns null then there was a numerical error calculating the category rates
                    // (probably a very small alpha) so reject the move.
                    return Double.NEGATIVE_INFINITY;
                }

                if (flip[k]) {
                    categoryRateBufferHelper[k].flipOffset(0);
                }

                beagle.setCategoryRatesWithIndex(categoryRateBufferHelper[k].getOffsetIndex(0), categoryRates);
                updatePartition[k] = true;
                if (DEBUG) {
                    System.out.println("updateSiteRateModels, updatePartition["+k+"] = " + updatePartition[k]);
                }
                updateAllPartitions = false;
            }
            k++;
        }

        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchLength();
            branchUpdateCount++;
        }

        if (branchUpdateCount > 0) {
            // TODO below only applies to homogenous substitution models

            int   [] eigenDecompositionIndices = new int   [branchUpdateCount * partitionCount];
            int   [] categoryRateIndices       = new int   [branchUpdateCount * partitionCount];
            int   [] probabilityIndices        = new int   [branchUpdateCount * partitionCount];
            double[] edgeLengths               = new double[branchUpdateCount * partitionCount];

            int op = 0;
            int partition = 0;
            for (EvolutionaryProcessDelegate evolutionaryProcessDelegate : evolutionaryProcessDelegates) {
                if (updatePartition[partition] || updateAllPartitions) {
                    if (flip[partition]) {
                        evolutionaryProcessDelegate.flipTransitionMatrices(branchUpdateIndices,
                                branchUpdateCount);
                    }

                    for (int i = 0; i < branchUpdateCount; i++) {
                        eigenDecompositionIndices[op] = evolutionaryProcessDelegate.getEigenIndex(0);
                        categoryRateIndices[op] = categoryRateBufferHelper[partition].getOffsetIndex(0);
                        probabilityIndices[op] = evolutionaryProcessDelegate.getMatrixIndex(branchUpdateIndices[i]);
                        edgeLengths[op] = branchLengths[i];
                        op++;
                    }
                }
                partition++;
            }

            beagle.updateTransitionMatricesWithMultipleModels(
                    eigenDecompositionIndices,
                    categoryRateIndices,
                    probabilityIndices,
                    null, // firstDerivativeIndices
                    null, // secondDerivativeIndices
                    edgeLengths,
                    op);
        }

        for (int i = 0; i < partitionCount; i++) {
            if (updatePartition[i] || updateAllPartitions) {

                if (DEBUG) {
                    System.out.println("updatePartition["+i+"] = " + updatePartition[i] + ", updateAllPartitions = " + updateAllPartitions);
                }

                if (flip[i]) {
                    // Flip all the buffers to be written to first...

                    for (NodeOperation op : nodeOperations) {
                        partialBufferHelper[i].flipOffset(op.getNodeNumber());
                    }
                }
            }
        }


        int operationCount = 0;
        k = 0;
        for (NodeOperation op : nodeOperations) {
            int nodeNum = op.getNodeNumber();

            int[] writeScale = new int[partitionCount];
            int[] readScale = new int[partitionCount];

            for (int i = 0; i < partitionCount; i++) {
                if (updatePartition[i] || updateAllPartitions) {
                    if (useScaleFactors[i]) {
                        // get the index of this scaling buffer
                        int n = nodeNum - tipCount;

                        if (recomputeScaleFactors[i]) {
                            // flip the indicator: can take either n or (internalNodeCount + 1) - n
                            scaleBufferHelper[i].flipOffset(n);

                            // store the index
                            scaleBufferIndices[i][n] = scaleBufferHelper[i].getOffsetIndex(n);

                            writeScale[i] = scaleBufferIndices[i][n]; // Write new scaleFactor
                            readScale[i] = Beagle.NONE;

                        } else {
                            writeScale[i] = Beagle.NONE;
                            readScale[i] = scaleBufferIndices[i][n]; // Read existing scaleFactor
                        }

                    } else {

                        writeScale[i] = Beagle.NONE; // Not using scaleFactors
                        readScale[i] = Beagle.NONE;
                    }
                }
            }

            //Example 1: 1 partition with 1 evolutionary model & -beagle_instances 3
            //partition 0 -> model 0
            //partition 1 -> model 0
            //partition 2 -> model 0

            //Example 2: 3 partitions with 3 evolutionary models & -beagle_instances 2
            //partitions 0 & 1 -> model 0
            //partitions 2 & 3 -> model 1
            //partitions 4 & 5 -> model 2

            int mapPartition = partitionCount / evolutionaryProcessDelegates.size();

            for (int i = 0; i < partitionCount; i++) {
                if (updatePartition[i] || updateAllPartitions) {

                    EvolutionaryProcessDelegate evolutionaryProcessDelegate = evolutionaryProcessDelegates.get(i / (mapPartition));
                    /*if (evolutionaryProcessDelegates.size() == partitionCount) {
                        evolutionaryProcessDelegate = evolutionaryProcessDelegates.get(i);
                    } else {
                        evolutionaryProcessDelegate = evolutionaryProcessDelegates.get(0);
                    }*/

                    operations[k] = partialBufferHelper[i].getOffsetIndex(nodeNum);
                    operations[k + 1] = writeScale[i];
                    operations[k + 2] = readScale[i];
                    operations[k + 3] = partialBufferHelper[i].getOffsetIndex(op.getLeftChild()); // source node 1
                    operations[k + 4] = evolutionaryProcessDelegate.getMatrixIndex(op.getLeftChild()); // source matrix 1
                    operations[k + 5] = partialBufferHelper[i].getOffsetIndex(op.getRightChild()); // source node 2
                    operations[k + 6] = evolutionaryProcessDelegate.getMatrixIndex(op.getRightChild()); // source matrix 2
                    operations[k + 7] = i;
                    //TODO: we don't know the cumulateScaleBufferIndex here yet (see below)
                    operations[k + 8] = Beagle.NONE;

                    if (DEBUG) {
                        if (k == 0 || (k == Beagle.PARTITION_OPERATION_TUPLE_SIZE)) {
                            System.out.println("write = " + writeScale[i] + "; read = " + readScale[i] + "; parent = " + operations[k] + ", k = " + k + ", i = " + i);
                        }

                    }

                    k += Beagle.PARTITION_OPERATION_TUPLE_SIZE;
                    operationCount++;
                }

            }
        }


        beagle.updatePartialsByPartition(operations, operationCount);

        //double[] rootPartials = new double[totalPatternCount * stateCount];
        //beagle.getPartials(rootIndex, 0, rootPartials);

        int[] cumulativeScaleIndices  = new int[partitionCount];
        for (int i = 0; i < partitionCount; i++) {
            cumulativeScaleIndices[i] = Beagle.NONE;
            if (useScaleFactors[i]) {
                if (recomputeScaleFactors[i] && (updatePartition[i] || updateAllPartitions)) {
                    scaleBufferHelper[i].flipOffset(internalNodeCount);
                    cumulativeScaleIndices[i] = scaleBufferHelper[i].getOffsetIndex(internalNodeCount);
                    //TODO: check with Daniel if calling these methods using an iteration can be done more efficiently
                    beagle.resetScaleFactorsByPartition(cumulativeScaleIndices[i], i);
                    beagle.accumulateScaleFactorsByPartition(scaleBufferIndices[i], internalNodeCount, cumulativeScaleIndices[i], i);
                } else {
                    cumulativeScaleIndices[i] = scaleBufferHelper[i].getOffsetIndex(internalNodeCount);
                }
            }
        }

//        double[] scaleFactors = new double[totalPatternCount];
//        beagle.getLogScaleFactors(cumulateScaleBufferIndex, scaleFactors);

        // these could be set only when they change but store/restore would need to be considered
        for (int i = 0; i < siteRateModels.size(); i++) {
            double[] categoryWeights = this.siteRateModels.get(i).getCategoryProportions();
            beagle.setCategoryWeights(i, categoryWeights);

            // This should probably explicitly be the state frequencies for the root node...
            double[] frequencies = evolutionaryProcessDelegates.get(i).getRootStateFrequencies();
            beagle.setStateFrequencies(i, frequencies);
        }


        if (DEBUG) {
            for (int i = 0; i < partitionCount; i++) {
                System.out.println("useScaleFactors=" + useScaleFactors[i] + " recomputeScaleFactors=" + recomputeScaleFactors[i] + " (partition: " + i + ")");
            }
        }

                /*System.out.println("partitionCount = " + partitionCount);
                for (int i = 0; i < partitionCount; i++) {
                    System.out.println("partitionIndices[" + i + "] = " + partitionIndices[i]);
                }*/

        int[] partitionIndices = new int[partitionCount];
        int[] rootIndices             = new int[partitionCount];
        int[] categoryWeightsIndices  = new int[partitionCount];
        int[] stateFrequenciesIndices = new int[partitionCount];
        double[] sumLogLikelihoods = new double[1];
        double[] sumLogLikelihoodsByPartition = new double[partitionCount];

        // create a list of partitions that have been updated
        int updatedPartitionCount = 0;
        for (int i = 0; i < partitionCount; i++) {
            if (updatePartition[i] || updateAllPartitions) {
                partitionIndices[updatedPartitionCount] = i;

                rootIndices            [updatedPartitionCount]  = partialBufferHelper[i].getOffsetIndex(rootNodeNumber);
                categoryWeightsIndices [updatedPartitionCount]  = i % siteRateModels.size();
                stateFrequenciesIndices[updatedPartitionCount]  = i % siteRateModels.size();
                cumulativeScaleIndices [updatedPartitionCount]  = cumulativeScaleIndices[i];

                updatedPartitionCount++;
            }
        }


        //TODO: check these arguments with Daniel
        //TODO: partitionIndices needs to be set according to which partitions need updating?
        beagle.calculateRootLogLikelihoodsByPartition(
                rootIndices,
                categoryWeightsIndices,
                stateFrequenciesIndices,
                cumulativeScaleIndices,
                partitionIndices,
                updatedPartitionCount,
                1,
                sumLogLikelihoodsByPartition,
                sumLogLikelihoods);

                /*System.out.println();
                for (int i = 0; i < partitionCount; i++) {
                    System.out.println("partition " + i + " lnL = " + sumLogLikelihoodsByPartition[i]);
                }*/

        // write updated partition likelihoods to the cached array
        for (int i = 0; i < updatedPartitionCount; i++) {
            cachedLogLikelihoodsByPartition[partitionIndices[i]] = sumLogLikelihoodsByPartition[i];

            // clear the global flags
            updatePartition[partitionIndices[i]] = false;
            recomputeScaleFactors[partitionIndices[i]] = false;

            partitionWasUpdated[partitionIndices[i]] = true;
        }

        double tmpLogL = sumLogLikelihoods[0];

        if (DEBUG) {
            for (int i = 0; i < partitionCount; i++) {
                System.out.println("partition " + i + ": " + cachedLogLikelihoodsByPartition[i] +
                        (partitionWasUpdated[i] ? " [updated]" : ""));
            }
        }

        // If these are needed...
//        if (patternLogLikelihoods == null) {
//            patternLogLikelihoods = new double[totalPatternCount];
//        }
//        beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        updateSubstitutionModels(false);
        updateSiteRateModels(false);
        updateAllPartitions = true;

        if (Double.isNaN(tmpLogL) || Double.isInfinite(tmpLogL)) {
            // one or more of the updated partitions has underflowed

            if (DEBUG) {
                System.out.println("Double.isNaN(logL) || Double.isInfinite(logL)");
            }

            for (int i = 0; i < updatedPartitionCount; i++) {
                if (Double.isNaN(sumLogLikelihoodsByPartition[i]) || Double.isInfinite(sumLogLikelihoodsByPartition[i])) {
                    everUnderflowed[partitionIndices[i]] = true;
                }
            }

            if (firstRescaleAttempt) {
                if (delayRescalingUntilUnderflow || rescalingScheme == PartialsRescalingScheme.DELAYED) {
                    // show a message but only every 1000 rescales
                    if (rescalingMessageCount % 1000 == 0) {
                        if (rescalingMessageCount > 0) {
                            Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood (" + rescalingMessageCount + " messages not shown).");
                        } else {
                            Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling... (" + getId() + ")");
                        }
                    }
                    rescalingMessageCount += 1;

                }

                for (int i = 0; i < updatedPartitionCount; i++) {
                    if (delayRescalingUntilUnderflow || rescalingScheme == PartialsRescalingScheme.DELAYED) {
                        if (Double.isNaN(sumLogLikelihoodsByPartition[i]) || Double.isInfinite(sumLogLikelihoodsByPartition[i])) {
                            useScaleFactors[partitionIndices[i]] = true;
                            recomputeScaleFactors[partitionIndices[i]] = true;
                            updatePartition[partitionIndices[i]] = true;

                            // turn off double buffer flipping so the next call overwrites the
                            // underflowed buffers. Flip will be turned on again in storeState for
                            // next step
                            flip[partitionIndices[i]] = false;
                            
                            updateAllPartitions = false;
                            if (DEBUG) {
                                System.out.println("Double.isNaN(logL) || Double.isInfinite(logL) (partition index: " + partitionIndices[i] + ")");
                            }
                        }

                    }
                }

                firstRescaleAttempt = false;

                throw new LikelihoodUnderflowException();
            }

            return Double.NEGATIVE_INFINITY;
        } else {

            for (int i = 0; i < updatedPartitionCount; i++) {
                if (partitionWasUpdated[partitionIndices[i]]) {
                    if (!this.delayRescalingUntilUnderflow || everUnderflowed[partitionIndices[i]]) {
                        if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                            if (!initialEvaluation) {
                                rescalingCount[partitionIndices[i]]++;
                            }
                        }
                    }
                    partitionWasUpdated[partitionIndices[i]] = false;
                }

                //TODO: probably better to only switch back those booleans that were actually altered
                recomputeScaleFactors[partitionIndices[i]] = false;
                flip[partitionIndices[i]] = true;
            }

            firstRescaleAttempt = true;
            initialEvaluation = false;
        }

        //********************************************************************
        double logL = 0.0;
        for (double l : cachedLogLikelihoodsByPartition) {
            logL += l;
        }

        return logL;
    }

    /*public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        // No need to rescale partials
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
    }*/

    /*private void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }*/

    @Override
    public void makeDirty() {
        updateSiteRateModels();
        updateSubstitutionModels();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof SiteRateModel) {
            updateSiteRateModel((SiteRateModel)model);
        } else if (model instanceof BranchModel) {
            updateSubstitutionModel((BranchModel)model);
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
        for (int i = 0; i < partitionCount; i++) {
            partialBufferHelper[i].storeState();
            categoryRateBufferHelper[i].storeState();
        }
        for (EvolutionaryProcessDelegate evolutionaryProcessDelegate : evolutionaryProcessDelegates) {
            evolutionaryProcessDelegate.storeState();
        }

        for (int i = 0; i < partitionCount; i++) {
            if (useScaleFactors[i] ) { // Only store when actually used
                scaleBufferHelper[i].storeState();
                System.arraycopy(scaleBufferIndices[i], 0, storedScaleBufferIndices[i], 0, scaleBufferIndices[i].length);
                //storedRescalingCount = rescalingCount;
            }

            // turn on double buffering flipping (may have been turned off to enable a rescale)
            flip[i] = true;
        }
        System.arraycopy(cachedLogLikelihoodsByPartition, 0, storedCachedLogLikelihoodsByPartition, 0, cachedLogLikelihoodsByPartition.length);
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restoreState() {

        for (int i = 0; i < partitionCount; i++) {
            partialBufferHelper[i].restoreState();
            categoryRateBufferHelper[i].restoreState();
        }
        for (EvolutionaryProcessDelegate evolutionaryProcessDelegate : evolutionaryProcessDelegates) {
            evolutionaryProcessDelegate.restoreState();
        }

        for (int i = 0; i < partitionCount; i++) {
            if (useScaleFactors[i]) {
                scaleBufferHelper[i].restoreState();
                int[] tmp = storedScaleBufferIndices[i];
                storedScaleBufferIndices[i] = scaleBufferIndices[i];
                scaleBufferIndices[i] = tmp;
                //rescalingCount = storedRescalingCount;
            }
        }

        double[] tmp = cachedLogLikelihoodsByPartition;
        cachedLogLikelihoodsByPartition = storedCachedLogLikelihoodsByPartition;
        storedCachedLogLikelihoodsByPartition = tmp;

    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
        // Callback not necessary
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void acceptState() {
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
        return "Using BEAGLE likelihood calculation library";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.AYRES_2012_BEAGLE);
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int nodeCount;
    private int tipCount;
    private int internalNodeCount;

    private int[] branchUpdateIndices;
    private double[] branchLengths;

    //provide per partition buffer indices
    private int[][] scaleBufferIndices;
    private int[][] storedScaleBufferIndices;

    private int[] operations;

    //allow flipping per partition
    private boolean[] flip;
    private BufferIndexHelper[] partialBufferHelper;
    private BufferIndexHelper[] scaleBufferHelper;
    private BufferIndexHelper[] categoryRateBufferHelper;

    private PartialsRescalingScheme rescalingScheme;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private boolean delayRescalingUntilUnderflow = true;

    //allow per partition rescaling
    private boolean[] useScaleFactors;
    private boolean[] recomputeScaleFactors;

    //keep track of underflow on a per partition basis
    private boolean[] everUnderflowed;

    private int[] rescalingCount;
    private int[] rescalingCountInner;

    private boolean firstRescaleAttempt;
    private int rescalingMessageCount = 0;

    // keep track of which partitions need to be updated
    private boolean updatePartition[];
    private boolean updateAllPartitions;
    private boolean partitionWasUpdated[];

    /**
     * the patternLists
     */
    private List<PatternList> patternLists;

    private final DataType dataType;

    private final int partitionCount;

    /**
     * the pattern weights across all patterns
     */
    private final double[] patternWeights;

    /**
     * The partition for each pattern
     */
    private final int[] patternPartitions;

    /**
     * the number of patterns for each partition
     */
    private final int[] patternCounts;

    /**
     * total number of patterns across all partitions
     */
    private final int totalPatternCount;

    /**
     * the number of states in the data
     */
    private final int stateCount;

    /**
     * the branch-site model for these sites
     */
    private final List<BranchModel> branchModels = new ArrayList<BranchModel>();

    /**
     * A delegate to handle substitution models on branches
     */
    private final List<EvolutionaryProcessDelegate> evolutionaryProcessDelegates = new ArrayList<EvolutionaryProcessDelegate>();

    /**
     * the site model for these sites
     */
    private final List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();

    /**
     * the number of rate categories
     */
    private final int categoryCount;

    /**
     * the BEAGLE library instance
     */
    private final Beagle beagle;

    /**
     * Cached log likelihood for each partition
     */
    private double[] cachedLogLikelihoodsByPartition;
    private double[] storedCachedLogLikelihoodsByPartition;

    /**
     * Flag to specify that the substitution model has changed
     */
    private final boolean[] updateSubstitutionModels;

    /**
     * Flag to specify that the site model has changed
     */
    private final boolean[] updateSiteRateModels;

    /**
     * Flag to take into account the first likelihood evaluation when initiating the MCMC chain
     */
    private boolean initialEvaluation = true;

}
