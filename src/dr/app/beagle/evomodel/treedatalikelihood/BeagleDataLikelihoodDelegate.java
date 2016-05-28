/*
 * BeagleDataLikelihoodDelegate.java
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
 * A DataLikelihoodDelegate that uses BEAGLE
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

import beagle.*;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.treelikelihood.*;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

public class BeagleDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate {
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
     *
     * @param tree Used for configuration - shouldn't be watched for changes
     * @param branchModel Specifies substitution model for each branch
     * @param patternList List of patterns
     * @param siteRateModel Specifies rates per site
     * @param useAmbiguities Whether to respect state ambiguities in data
     * @param rescalingScheme BEAGLE rescaling scheme
     * @param delayRescalingUntilUnderflow whether to delay rescaling until first underflow
     */
    public BeagleDataLikelihoodDelegate(Tree tree,
                                        BranchModel branchModel,
                                        PatternList patternList,
                                        SiteRateModel siteRateModel,
                                        boolean useAmbiguities,
                                        PartialsRescalingScheme rescalingScheme,
                                        boolean delayRescalingUntilUnderflow) {

        super("BeagleDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("Using BEAGLE DataLikelihood Delegate");

        this.patternList = patternList;
        this.dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

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

        try {


            int compactPartialsCount = tipCount;
            if (useAmbiguities) {
                // if we are using ambiguities then we don't use tip partials
                compactPartialsCount = 0;
            }

            // one partials buffer for each tip and two for each internal node (for store restore)
            partialBufferHelper = new dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper(nodeCount, tipCount);

            // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
            scaleBufferHelper = new dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper(getScaleBufferCount(), 0);

            substitutionModelDelegate = new SubstitutionModelDelegate(tree, branchModel, -1);

            // first set the rescaling scheme to use from the parser
            this.rescalingScheme = rescalingScheme;
            this.delayRescalingUntilUnderflow = delayRescalingUntilUnderflow;

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


            if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
                if (stateCount == 4 && patternList.getPatternCount() < 10000) // TODO determine good cut-off
                    preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
            }

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

            if (substitutionModelDelegate.canReturnComplexDiagonalization()) {
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

            if (patternList instanceof UncertainSiteList) {
                useAmbiguities = true;
            }

            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");

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

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
        this.useAmbiguities = useAmbiguities;
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
                                   PatternList patternList,
                                   int sequenceIndex,
                                   int nodeIndex) {
        double[] partials = new double[patternCount * stateCount * categoryCount];

        boolean[] stateSet;

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

            if (patternList instanceof UncertainSiteList) {
                ((UncertainSiteList) patternList).fillPartials(sequenceIndex, i, partials, v);
                v += stateCount;
                // TODO Add this functionality to SimpleSiteList to avoid if statement here
            } else {

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


    //    public void setStates(int tipIndex, int[] states) {
//        System.err.println("BTL:setStates");
//        beagle.setTipStates(tipIndex, states);
//        makeDirty();
//    }
//
//    public void getStates(int tipIndex, int[] states) {
//        System.err.println("BTL:getStates");
//        beagle.getTipStates(tipIndex, states);
//    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {

        recomputeScaleFactors = false;

        if (!this.delayRescalingUntilUnderflow || everUnderflowed) {
            if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS || this.rescalingScheme == PartialsRescalingScheme.DELAYED) {
                useScaleFactors = true;
                recomputeScaleFactors = true;
            } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                useScaleFactors = true;

                if (rescalingCount > rescalingFrequency) {
                    rescalingCount = 0;
                    rescalingCountInner = 0;
                }

                if (rescalingCountInner < RESCALE_TIMES) {
                    recomputeScaleFactors = true;
                    updateAllNodes();

                    rescalingCountInner++;
                }

                rescalingCount++;
            }
        }

        branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
           branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchNumber();
            branchUpdateCount ++;
        }

        if (updateSubstitutionModel) { // TODO More efficient to update only the substitution model that changed, instead of all
            substitutionModelDelegate.updateSubstitutionModels(beagle);

            // we are currently assuming a no-category model...
        }

        if (updateSiteModel) {
            double[] categoryRates = this.siteRateModel.getCategoryRates();
            beagle.setCategoryRates(categoryRates);
        }

        if (branchUpdateCount > 0) {
            substitutionModelDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount);
        }

        operationCount = 0;
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            operations[k] = op.getNodeNumber();
            operations[k + 1] = op.getLeftChild();
            operations[k + 2] = op.getRightChild();
            operationCount++;
        }


        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

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

            logL = sumLogLikelihoods[0];

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (delayRescalingUntilUnderflow || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling
                    if (rescalingScheme == PartialsRescalingScheme.DYNAMIC || (rescalingCount == 0)) {
                        Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling...");
                    }
                    useScaleFactors = true;
                    recomputeScaleFactors = true;

                    branchUpdateCount = 0;

                    updateAllNodes();

                    if (hasRestrictedPartials) {
                        for (int i = 0; i <= numRestrictedPartials; i++) {
                            operationCount[i] = 0;
                        }
                    } else {
                        operationCount[0] = 0;
                    }

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(treeModel, root, null, false);

                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

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
        updateSiteModel = true;
        updateSubstitutionModel = true;
        updateAllPatterns();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    /**
     * Set update flag for a pattern
     */
    private void updatePattern(int i) {
        if (updatePattern != null) {
            updatePattern[i] = true;
        }
    }

    /**
     * Set update flag for all patterns
     */
    private void updateAllPatterns() {
        if (updatePattern != null) {
            for (int i = 0; i < patternCount; i++) {
                updatePattern[i] = true;
            }
        }
    }

    /**
     * Stores the additional state other than model components
     */
    public void storeState() {
        partialBufferHelper.storeState();
        substitutionModelDelegate.storeState();

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }
    }

    /**
     * Restore the additional stored state
     */
    public void restoreState() {
        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore

        partialBufferHelper.restoreState();
        substitutionModelDelegate.restoreState();

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
    private int branchUpdateCount;

    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[] operations;
    private int operationCount;

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
     * the patternList
     */
    private PatternList patternList = null;
    private DataType dataType = null;

    /**
     * the pattern weights
     */
    private double[] patternWeights;

    /**
     * the number of patterns
     */
    private int patternCount;

    /**
     * the number of states in the data
     */
    private int stateCount;

    /**
     * Flags to specify which patterns are to be updated
     */
    private boolean[] updatePattern = null;

    /**
     * the branch-site model for these sites
     */
    private final BranchModel branchModel;

    /**
     * A delegate to handle substitution models on branches
     */
    private final SubstitutionModelDelegate substitutionModelDelegate;

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
    private boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    private boolean updateSiteModel;

//    /***
//     * Flag to specify if LikelihoodCore supports dynamic rescaling
//     */
//    private boolean dynamicRescaling = false;


    /**
     * Flag to specify if site patterns are acertained
     */

    private boolean ascertainedSitePatterns = false;


    /**
     * Flag to specify if ambiguity codes are in use
     */
    private final boolean useAmbiguities;
}
