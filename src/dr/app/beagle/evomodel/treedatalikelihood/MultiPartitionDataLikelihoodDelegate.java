/*
 * MultiPartitionDataLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.evomodel.treedatalikelihood;/**
 * BeagleDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate that uses BEAGLE 3 to allow for parallelization across multiple data partitions
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

import beagle.*;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.EpochBranchModel;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.treelikelihood.SubstitutionModelDelegate;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

public class MultiPartitionDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate {
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
                                                boolean useAmbiguities) {

        super("MultiPartitionDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("Using Multi-Partition Data Likelihood Delegate");

        this.patternLists.addAll(patternLists);
        this.dataType = this.patternLists.get(0).getDataType();
        patternCounts = new int[this.patternLists.size()];
        patternWeights = new double[this.patternLists.size()][];
        int k = 0;
        for (PatternList patternList : this.patternLists) {
            assert(patternList.getDataType().equals(this.dataType));
            patternCounts[k] = patternList.getPatternCount();
            patternWeights[k] = patternList.getPatternWeights();
            k++;
        }

        stateCount = dataType.getStateCount();

        this.useAmbiguities = useAmbiguities;

        assert(branchModels.size() == 1 || branchModels.size() == patternLists.size());
        this.branchModels.addAll(branchModels);
        for (BranchModel branchModel : this.branchModels) {
            addModel(branchModel);
        }

        assert(siteRateModels.size() == 1 || siteRateModels.size() == patternLists.size());
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
        scaleBufferIndices = new int[internalNodeCount];
        storedScaleBufferIndices = new int[internalNodeCount];

        operations = new int[internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];

        try {


            int compactPartialsCount = tipCount;
            if (useAmbiguities) {
                // if we are using ambiguities then we don't use tip partials
                compactPartialsCount = 0;
            }

            // one partials buffer for each tip and two for each internal node (for store restore)
            partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

            // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
            scaleBufferHelper = new BufferIndexHelper(getScaleBufferCount(), 0);

            // create a substitutionModelDelegate for each branchModel
            for (BranchModel branchModel : this.branchModels) {
                substitutionModelDelegates.add(new SubstitutionModelDelegate(tree, branchModel, -1));
            }

            // first set the rescaling scheme to use from the parser
            this.rescalingScheme = PartialsRescalingScheme.ALWAYS;
            this.delayRescalingUntilUnderflow = false;

            int[] resourceList = null;
            long preferenceFlags = 0;
            long requirementFlags = 0;

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

            if (substitutionModelDelegates.get(0).canReturnComplexDiagonalization()) {
                requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
            }

            beagle = BeagleFactory.loadBeagleInstance(
                    tipCount,
                    partialBufferHelper.getBufferCount(),
                    compactPartialsCount,
                    stateCount,
                    patternCount,
                    substitutionModelDelegate.getEigenBufferCount(),
                    substitutionModelDelegate.getMatrixBufferCount(),
                    categoryCount,
                    scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
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

            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            String patternCountString = "" + this.patternLists.get(0).getPatternCount();
            for (int i = 1; i < this.patternLists.size(); i++) {
                patternCountString += ", " + patternLists.get(i).getPatternCount();
            }
            logger.info("  With " + patternLists.size() + " partitions comprising " + patternCountString + " unique site patterns");

            k = 0;
            for (PatternList patternList : this.patternLists) {
                // @todo - should check that each patternList spans the same set of taxa
                for (int i = 0; i < tipCount; i++) {
                    // Find the id of tip i in the patternList
                    String id = tree.getTaxonId(i);
                    int index = patternList.getTaxonIndex(id);

                    if (index == -1) {
                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + tree.getId() +
                                ", is not found in patternList, " + patternList.getId());
                    } else {
                        if (useAmbiguities) {
                            setPartials(beagle, k, patternList, index, i);
                        } else {
                            setStates(beagle, k, patternList, index, i);
                        }
                    }
                }

                beagle.setPatternWeights(k, patternWeights[k]);

                k++;
            }

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

            updateSubstitutionModels = new boolean[substitutionModelDelegates.size()];
            updateSubstitutionModels();

            updateSiteModels = new boolean[branchModels.size()];
            updateSiteModels();

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
    }

    private void updateSubstitutionModels() {
        for (int i = 0; i < substitutionModelDelegates.size(); i++) {
            updateSubstitutionModels[i] = true;
        }
    }

    private void updateSubstitutionModel(SubstitutionModelDelegate substitutionModelDelegate) {
        for (int i = 0; i < substitutionModelDelegates.size(); i++) {
            if (substitutionModelDelegates.get(i) == substitutionModelDelegate) {
                updateSubstitutionModels[i] = true;
            }
        }
    }

    private void updateSiteModels() {
        for (int i = 0; i < branchModels.size(); i++) {
            updateSiteModels[i] = true;
        }
    }

    private void updateSiteModel(BranchModel branchModel) {
        for (int i = 0; i < branchModels.size(); i++) {
            if (branchModels.get(i) == branchModel) {
                updateSiteModels[i] = true;
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
     * @param patternList   patternList
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    private final void setPartials(Beagle beagle,
                                   int partitionIndex,
                                   PatternList patternList,
                                   int sequenceIndex,
                                   int nodeIndex) {

        int patternCount = patternList.getPatternCount();
        double[] partials = new double[patternCount * stateCount * categoryCount];

        boolean[] stateSet;

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

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

        // if there is more than one category then replicate the partials for each
        int n = patternCount * stateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(partitionIndex, nodeIndex, partials);
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
                                 int partitionIndex,
                                 PatternList patternList,
                                 int sequenceIndex,
                                 int nodeIndex) {

        int patternCount = patternList.getPatternCount();

        int i;

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {

            states[i] = patternList.getPatternState(sequenceIndex, i);
        }

        beagle.setTipStates(partitionIndex, nodeIndex, states);
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) throws LikelihoodUnderflowException {

        // For the first version just do scaling always
        useScaleFactors = true;
        recomputeScaleFactors = true;

        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchNumber();
            branchUpdateCount++;
        }

        int k = 0;
        for (SubstitutionModelDelegate substitutionModelDelegate : substitutionModelDelegates) {
            if (updateSubstitutionModels[k]) {
                // TODO More efficient to update only the substitution model that changed, instead of all
                substitutionModelDelegate.updateSubstitutionModels(beagle);

                // we are currently assuming a no-category model...
            }
            k++;
        }

         k = 0;
        for (SiteRateModel siteRateModel : siteRateModels) {
            if (updateSiteModels[k]) {
                double[] categoryRates = siteRateModelgetCategoryRates();
                beagle.setCategoryRates(k, categoryRates);
            }
            k++;
        }

        if (branchUpdateCount > 0) {
            substitutionModelDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount);
        }

        int operationCount = nodeOperations.size();
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            int nodeNum = op.getNodeNumber();

            if (flip) {
                // first flip the partialBufferHelper
                partialBufferHelper.flipOffset(nodeNum);
            }

            operations[k] = partialBufferHelper.getOffsetIndex(nodeNum);

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
            operations[k + 4] = substitutionModelDelegate.getMatrixIndex(op.getLeftChild()); // source matrix 1
            operations[k + 5] = partialBufferHelper.getOffsetIndex(op.getRightChild()); // source node 2
            operations[k + 6] = substitutionModelDelegate.getMatrixIndex(op.getRightChild()); // source matrix 2

            k += Beagle.OPERATION_TUPLE_SIZE;
        }

        beagle.updatePartials(operations, operationCount, Beagle.NONE);

        int rootIndex = partialBufferHelper.getOffsetIndex(rootNodeNumber);

        double[] categoryWeights = this.siteRateModel.getCategoryProportions();

        // This should probably explicitly be the state frequencies for the root node...
        double[] frequencies = substitutionModelDelegate.getRootStateFrequencies();

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

        // these could be set only when they change but store/restore would need to be considered
        beagle.setCategoryWeights(0, categoryWeights);
        beagle.setStateFrequencies(0, frequencies);

        double[] sumLogLikelihoods = new double[1];

        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

        double logL = sumLogLikelihoods[0];

        if (Double.isNaN(logL) || Double.isInfinite(logL)) {
            everUnderflowed = true;
            // turn off double buffer flipping so the next call overwrites the
            // underflowed buffers. Flip will be turned on again in storeState for
            // next step
            flip = false;
            throw new LikelihoodUnderflowException();
        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
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
        updateSiteModels();
        updateSubstitutionModels();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof SiteRateModel) {
            updateSiteModel((SiteRateModel)model);
        } else if (model instanceof BranchModel) {
            updateSubstitutionModel((BranchModel)model);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    /**
     * Stores the additional state other than model components
     */
    public void storeState() {
        partialBufferHelper.storeState();
        for (SubstitutionModelDelegate substitutionModelDelegate : substitutionModelDelegates) {
            substitutionModelDelegate.storeState();
        }

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }

        // turn on double buffering flipping (may have been turned off to enable a rescale)
        flip = true;
    }

    /**
     * Restore the additional stored state
     */
    public void restoreState() {
        updateSiteModels(); // this is required to upload the categoryRates to BEAGLE after the restore

        partialBufferHelper.restoreState();
        for (SubstitutionModelDelegate substitutionModelDelegate : substitutionModelDelegates) {
            substitutionModelDelegate.restoreState();
        }

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp;
//            rescalingCount = storedRescalingCount;
        }

    }

    @Override
    protected void acceptState() {
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int nodeCount;
    private int tipCount;
    private int internalNodeCount;

    private int[] branchUpdateIndices;
    private double[] branchLengths;

    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[] operations;

    private boolean flip = true;
    private BufferIndexHelper partialBufferHelper;
    private BufferIndexHelper scaleBufferHelper;

    private PartialsRescalingScheme rescalingScheme;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private boolean delayRescalingUntilUnderflow = true;

    private boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    /**
     * the patternLists
     */
    private List<PatternList> patternLists = new ArrayList<PatternList>();
    private DataType dataType = null;

    /**
     * the pattern weights
     */
    private double[][] patternWeights;

    /**
     * the number of patterns
     */
    private int[] patternCounts;

    /**
     * the number of states in the data
     */
    private int stateCount;

    /**
     * the branch-site model for these sites
     */
    private final List<BranchModel> branchModels = new ArrayList<BranchModel>();

    /**
     * A delegate to handle substitution models on branches
     */
    private final List<SubstitutionModelDelegate> substitutionModelDelegates = new ArrayList<SubstitutionModelDelegate>();

    /**
     * the site model for these sites
     */
    private final List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();

    /**
     * the pattern likelihoods
     */
    private double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    private int categoryCount;

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
    private Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    private boolean[] updateSubstitutionModels;

    /**
     * Flag to specify that the site model has changed
     */
    private boolean[] updateSiteModels;

    public static void main(String[] args) {

        try {

            MathUtils.setSeed(666);

            System.out.println("Test case 1: simulateOnePartition");

            int sequenceLength = 1000;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            Tree tree = importer.importTree(null);
            TreeModel treeModel = new TreeModel(tree);

            // create Frequency Model
            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25,
                    0.25, 0.25});
            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                    freqs);

            // create branch model
            Parameter kappa1 = new Parameter.Default(1, 1);
            Parameter kappa2 = new Parameter.Default(1, 1);

            HKY hky1 = new HKY(kappa1, freqModel);
            HKY hky2 = new HKY(kappa2, freqModel);

            HomogeneousBranchModel homogenousBranchSubstitutionModel = new HomogeneousBranchModel(hky1);

            List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
            substitutionModels.add(hky1);
            substitutionModels.add(hky2);
            List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>();
            freqModels.add(freqModel);

            Parameter epochTimes = new Parameter.Default(1, 20);

            // create branch rate model
            Parameter rate = new Parameter.Default(1, 0.001);
            BranchRateModel branchRateModel = new StrictClockBranchRates(rate);

            // create site model
            List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
            siteRateModels.add(siteRateModel);

            List<BranchModel> branchModels = new ArrayList<BranchModel>();
            branchModels.add(new HomogeneousBranchModel(hky1));

//            BranchModel epochBranchModel = new EpochBranchModel(treeModel, substitutionModels, epochTimes);

            // create partition
            Partition partition1 = new Partition(treeModel, //
                    homogenousBranchSubstitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    sequenceLength - 1, // to
                    1 // every
            );

            partitionsList.add(partition1);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(partitionsList);
            Alignment alignment = simulator.simulate(false, false);

            List<PatternList> partitions = new ArrayList<PatternList>();
            partitions.add(alignment);

            MultiPartitionDataLikelihoodDelegate bdld = new MultiPartitionDataLikelihoodDelegate(treeModel, partitions, branchModels, siteRateModels, false);
            TreeDataLikelihood tdl  = new TreeDataLikelihood(bdld, treeModel, branchRateModel);

            System.out.println("BeagleDataLikelihoodDelegate(homogeneous) = " + tdl.getLogLikelihood());

//            nbtl = new BeagleTreeLikelihood(alignment, treeModel, epochBranchModel, siteRateModel, branchRateModel, null, false, PartialsRescalingScheme.DEFAULT, false);
//
//            System.out.println("nBTL(epoch) = " + nbtl.getLogLikelihood());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block
    }

}
