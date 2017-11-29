/*
 * BeagleTreeLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treelikelihood;

import beagle.*;
import dr.evolution.datatype.HiddenDataType;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.substmodel.MarkovModulatedSubstitutionModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodelxml.treelikelihood.BeagleTreeLikelihoodParser;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ThreadAwareLikelihood;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;
import java.util.logging.Logger;

/**
 * BeagleTreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 * @version $Id$
 */

@SuppressWarnings("serial")
@Deprecated // Switching to TreeDataLikelihood
public class BeagleTreeLikelihood extends AbstractSinglePartitionTreeLikelihood implements ThreadAwareLikelihood, Citable {
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

    private static final boolean RESCALING_OFF = false; // a debugging switch
    private static final boolean DEBUG = false;

    public BeagleTreeLikelihood(PatternList patternList,
                                TreeModel treeModel,
                                BranchModel branchModel,
                                SiteRateModel siteRateModel,
                                BranchRateModel branchRateModel,
                                TipStatesModel tipStatesModel,
                                boolean useAmbiguities,
                                PartialsRescalingScheme rescalingScheme,
                                boolean delayRescalingUntilUnderflow) {

        this(patternList, treeModel, branchModel, siteRateModel, branchRateModel, tipStatesModel, useAmbiguities, rescalingScheme, delayRescalingUntilUnderflow, null);
    }

    public BeagleTreeLikelihood(PatternList patternList,
                                TreeModel treeModel,
                                BranchModel branchModel,
                                SiteRateModel siteRateModel,
                                BranchRateModel branchRateModel,
                                TipStatesModel tipStatesModel,
                                boolean useAmbiguities,
                                PartialsRescalingScheme rescalingScheme,
                                boolean delayRescalingUntilUnderflow,
                                Map<Set<String>, Parameter> partialsRestrictions) {

        super(BeagleTreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);

        try {
            final Logger logger = Logger.getLogger("dr.evomodel");

            logger.info("\nUsing BEAGLE TreeLikelihood");

            this.siteRateModel = siteRateModel;
            addModel(this.siteRateModel);

            this.branchModel = branchModel;
            addModel(this.branchModel);

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                logger.info("  Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultBranchRateModel();
            }
            addModel(this.branchRateModel);

            if (patternList instanceof UncertainSiteList ||
                    patternList.getDataType() instanceof HiddenDataType) {
                useAmbiguities = true;
            }

            this.tipStatesModel = tipStatesModel;

            this.categoryCount = this.siteRateModel.getCategoryCount();

            this.tipCount = treeModel.getExternalNodeCount();

            internalNodeCount = nodeCount - tipCount;

            int compactPartialsCount = tipCount;
            if (useAmbiguities) {
                // if we are using ambiguities then we don't use tip partials
                compactPartialsCount = 0;
            }

            // one partials buffer for each tip and two for each internal node (for store restore)
            partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

            // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
            scaleBufferHelper = new BufferIndexHelper(getScaleBufferCount(), 0);

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

            int extraBufferCount = -1; // default
            if (extraBufferOrder.size() > 0) {
                extraBufferCount = extraBufferOrder.get(instanceCount % extraBufferOrder.size());
            }
            substitutionModelDelegate = new SubstitutionModelDelegate(treeModel, branchModel, extraBufferCount);

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

            // TODO Remove once issue #854 is fixed
            for (int s = 0; s < substitutionModelDelegate.getSubstitutionModelCount(); ++s) {
                if (substitutionModelDelegate.getSubstitutionModel(s) instanceof MarkovModulatedSubstitutionModel) {
                    this.rescalingScheme = PartialsRescalingScheme.ALWAYS;
                    this.delayRescalingUntilUnderflow = false;
                }
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

            // Check for matching state counts
            int stateCount2 = branchModel.getRootFrequencyModel().getFrequencyCount();
            if (stateCount != stateCount2) {
                throw new RuntimeException("Pattern state count (" + stateCount
                        + ") does not match substitution model state count (" + stateCount2 + ")");
            }

            instanceCount++;

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
            logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");

            if (tipStatesModel != null) {
                tipStatesModel.setTree(treeModel);

                if (tipStatesModel.getModelType() == TipStatesModel.Type.PARTIALS) {
                    tipPartials = new double[patternCount * stateCount];
                } else {
                    tipStates = new int[patternCount];
                }

                addModel(tipStatesModel);
            }

            for (int i = 0; i < tipCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                } else {
                    if (tipStatesModel != null) {
                        // using a tipPartials model.
                        // First set the observed states:
                        tipStatesModel.setStates(patternList, index, i, id);

                        if (tipStatesModel.getModelType() == TipStatesModel.Type.PARTIALS) {
                            // Then set the tip partials as determined by the model:
                            setPartials(beagle, tipStatesModel, i);
                        } else {
                            // or the tip states:
                            tipStatesModel.getTipStates(i, tipStates);
                            beagle.setTipStates(i, tipStates);
                        }

                    } else {
                        if (useAmbiguities) {
                            setPartials(beagle, patternList, index, i);
                        } else {
                            setStates(beagle, patternList, index, i);
                        }
                    }
                }
            }

            if (patternList instanceof AscertainedSitePatterns) {
                ascertainedSitePatterns = true;
            }

            this.partialsRestrictions = partialsRestrictions;
//            hasRestrictedPartials = (partialsRestrictions != null);
            if (hasRestrictedPartials) {
                numRestrictedPartials = partialsRestrictions.size();
                updateRestrictedNodePartials = true;
                partialsMap = new Parameter[treeModel.getNodeCount()];
                partials = new double[stateCount * patternCount * categoryCount];
            } else {
                numRestrictedPartials = 0;
                updateRestrictedNodePartials = false;
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
                rescaleMessage += (parenthesis ? ", " : " (") + "delay rescaling until first overflow";
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
        hasInitialized = true;
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

    public TipStatesModel getTipStatesModel() {
        return tipStatesModel;
    }

    public PatternList getPatternsList() {
        return patternList;
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public BranchModel getBranchModel() {
        return branchModel;
    }

    public SiteRateModel getSiteRateModel() {
        return siteRateModel;
    }

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    public PartialsRescalingScheme getRescalingScheme() {
        return rescalingScheme;
    }

    public boolean isDelayRescalingUntilUnderflow() {
        return delayRescalingUntilUnderflow;
    }

    public Map<Set<String>, Parameter> getPartialsRestrictions() {
        return partialsRestrictions;
    }

    public boolean useAmbiguities() {
        return useAmbiguities;
    }

    protected int getScaleBufferCount() {
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
    protected final void setPartials(Beagle beagle,
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
    protected final void setPartials(Beagle beagle,
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

    public int getPatternCount() {
        return patternCount;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternList   patternList
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    protected final void setStates(Beagle beagle,
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

    public final void setPatternWeights(double[] patternWeights) {
        this.patternWeights = patternWeights;
        beagle.setPatternWeights(patternWeights);
    }


    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        fireModelChanged();

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeChangedEvent) object).getNode());
                    updateRestrictedNodePartials = true;

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                    updateRestrictedNodePartials = true;
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else if (model == branchModel) {
            if (index == -1) {
                updateSubstitutionModel = true;
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }
        } else if (model == siteRateModel) {

            updateSiteModel = true;
            updateAllNodes();

        } else if (model == tipStatesModel) {
            if (object instanceof Taxon) {
                for (int i = 0; i < treeModel.getNodeCount(); i++)
                    if (treeModel.getNodeTaxon(treeModel.getNode(i)) != null && treeModel.getNodeTaxon(treeModel.getNode(i)).getId().equalsIgnoreCase(((Taxon) object).getId()))
                        updateNode(treeModel.getNode(i));
            } else if (object instanceof Parameter) {
                // ignore...
            } else {
                updateAllNodes();
            }
        } else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }

        super.handleModelChangedEvent(model, object, index);
    }

    @Override
    public void makeDirty() {
        super.makeDirty();
        updateSiteModel = true;
        updateSubstitutionModel = true;
        updateRestrictedNodePartials = true;
    }
// **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {
        partialBufferHelper.storeState();
        substitutionModelDelegate.storeState();

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }

        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {
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

        updateRestrictedNodePartials = true;

        super.restoreState();

    }

//    int marcCount = 0;
    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (branchUpdateIndices == null) {
            branchUpdateIndices = new int[nodeCount];
            branchLengths = new double[nodeCount];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[numRestrictedPartials + 1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[numRestrictedPartials + 1];
        }

        recomputeScaleFactors = false;

        if (!this.delayRescalingUntilUnderflow || everUnderflowed) {
            if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS || this.rescalingScheme == PartialsRescalingScheme.DELAYED) {
                useScaleFactors = true;
                recomputeScaleFactors = true;
            } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
                useScaleFactors = true;
                if (DEBUG) {
                    System.out.println("rescalingCount = " + rescalingCount);
                    System.out.println("rescalingCountInner = " + rescalingCountInner);
                }

                if (rescalingCount > rescalingFrequency) {
                    if (DEBUG) {
                        System.out.println("rescalingCount > rescalingFrequency");
                    }
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

        if (RESCALING_OFF) { // a debugging switch
            useScaleFactors = false;
            recomputeScaleFactors = false;
        }

        if (tipStatesModel != null) {
            int tipCount = treeModel.getExternalNodeCount();
            for (int index = 0; index < tipCount; index++) {
                if (updateNode[index]) {
                    if (tipStatesModel.getModelType() == TipStatesModel.Type.PARTIALS) {
                        tipStatesModel.getTipPartials(index, tipPartials);
                        beagle.setTipPartials(index, tipPartials);
                    } else {
                        tipStatesModel.getTipStates(index, tipStates);
                        beagle.setTipStates(index, tipStates);
                    }
                }
            }
        }

        branchUpdateCount = 0;
        operationListCount = 0;

        if (hasRestrictedPartials) {
            for (int i = 0; i <= numRestrictedPartials; i++) {
                operationCount[i] = 0;
            }
        } else {
            operationCount[0] = 0;
        }

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root, null, true);

        if (DEBUG) {
            System.out.println("operationCount = " + operationCount[operationListCount]);
        }

        if (updateSubstitutionModel) { // TODO More efficient to update only the substitution model that changed, instead of all
            substitutionModelDelegate.updateSubstitutionModels(beagle);

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
        }

        if (branchUpdateCount > 0) {
            substitutionModelDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount);
        }

        if (COUNT_TOTAL_OPERATIONS) {
            totalMatrixUpdateCount += branchUpdateCount;

            for (int i = 0; i <= numRestrictedPartials; i++) {
                totalOperationCount += operationCount[i];
            }
        }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            if (hasRestrictedPartials) {
                for (int i = 0; i <= numRestrictedPartials; i++) {
                    beagle.updatePartials(operations[i], operationCount[i], Beagle.NONE);
                    if (i < numRestrictedPartials) {
//                        restrictNodePartials(restrictedIndices[i]);
                    }
                }
            } else {
                beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);
            }

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNumber());

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

            if (DEBUG) {
                System.out.println("useScaleFactors=" + useScaleFactors + " recomputeScaleFactors=" + recomputeScaleFactors);
            }

            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            logL = sumLogLikelihoods[0];

            if (DEBUG) {
                System.out.println(logL);
                // if (logL > -90000) {
                //     System.exit(0);
                // }
            }

            beagle.getSiteLogLikelihoods(patternLogLikelihoods);

            if (ascertainedSitePatterns) {
                // Need to correct for ascertainedSitePatterns
                beagle.getSiteLogLikelihoods(patternLogLikelihoods);
                logL = getAscertainmentCorrectedLogLikelihood((AscertainedSitePatterns) patternList,
                        patternLogLikelihoods, patternWeights);
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {

                if (DEBUG) {
                    System.out.println("Double.isNaN(logL) || Double.isInfinite(logL)");
                }

                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (delayRescalingUntilUnderflow || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling
                    if (rescalingScheme == PartialsRescalingScheme.DYNAMIC || (rescalingCount == 0)) {
                        // show a message but only every 1000 rescales
                        if (rescalingMessageCount % 1000 == 0) {
                            if (rescalingMessageCount > 0) {
                                Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood (" + rescalingMessageCount + " messages not shown).");
                            } else {
                                Logger.getLogger("dr.evomodel").info("Underflow calculating likelihood. Attempting a rescaling...");
                            }
                        }
                        rescalingMessageCount += 1;
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

        // If these are needed...
        //beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
        //********************************************************************

        return logL;
    }

    public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        /* No need to rescale partials */
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number), cumulativeBufferIndex, partials);
    }

    public boolean arePartialsRescaled() {
        return useScaleFactors;
    }

    protected void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }

    private void restrictNodePartials(int nodeIndex) {

        Parameter restrictionParameter = partialsMap[nodeIndex];
        if (restrictionParameter == null) {
            return;
        }

        getPartials(nodeIndex, partials);

        double[] restriction = restrictionParameter.getParameterValues();
        final int partialsLengthPerCategory = stateCount * patternCount;
        if (restriction.length == partialsLengthPerCategory) {
            for (int i = 0; i < categoryCount; i++) {
                componentwiseMultiply(partials, partialsLengthPerCategory * i, restriction, 0, partialsLengthPerCategory);
            }
        } else {
            componentwiseMultiply(partials, 0, restriction, 0, partialsLengthPerCategory * categoryCount);
        }

        setPartials(nodeIndex, partials);
    }

    private void componentwiseMultiply(double[] a, final int offsetA, double[] b, final int offsetB, final int length) {
        for (int i = 0; i < length; i++) {
            a[offsetA + i] *= b[offsetB + i];
        }
    }

    private void computeNodeToRestrictionMap() {
        Arrays.fill(partialsMap, null);
        for (Set<String> taxonNames : partialsRestrictions.keySet()) {
            NodeRef node = TreeUtils.getCommonAncestorNode(treeModel, taxonNames);
            partialsMap[node.getNumber()] = partialsRestrictions.get(taxonNames);
        }
    }

    private double getAscertainmentCorrectedLogLikelihood(AscertainedSitePatterns patternList,
                                                          double[] patternLogLikelihoods,
                                                          double[] patternWeights) {
        double logL = 0.0;
        double ascertainmentCorrection = patternList.getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }
        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param tree           tree
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private boolean traverse(Tree tree, NodeRef node, int[] operatorNumber, boolean flip) {

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate;

            synchronized (branchRateModel) {
                branchRate = branchRateModel.getBranchRate(tree, node);
            }
            final double parentHeight = tree.getNodeHeight(parent);
            final double nodeHeight = tree.getNodeHeight(node);

            // Get the operational time of the branch
            final double branchLength = branchRate * (parentHeight - nodeHeight);
            if (branchLength < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchLength + " (parent: " + parent +
                        "; height: " + parentHeight + " - child: " + node + "height: " + nodeHeight + ")");
            }

            if (flip) {
                substitutionModelDelegate.flipMatrixBuffer(nodeNum);
            }
            branchUpdateIndices[branchUpdateCount] = nodeNum;
            branchLengths[branchUpdateCount] = branchLength;
            branchUpdateCount++;

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final int[] op1 = {-1};
            final boolean update1 = traverse(tree, child1, op1, flip);

            NodeRef child2 = tree.getChild(node, 1);
            final int[] op2 = {-1};
            final boolean update2 = traverse(tree, child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;

                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);

                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;

                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }

                } else {

                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNumber()); // source node 1
                operations[x + 4] = substitutionModelDelegate.getMatrixIndex(child1.getNumber()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNumber()); // source node 2
                operations[x + 6] = substitutionModelDelegate.getMatrixIndex(child2.getNumber()); // source matrix 2

                operationCount[operationListCount]++;

                update = true;

                if (hasRestrictedPartials) {
                    // Test if this set of partials should be restricted
                    if (updateRestrictedNodePartials) {
                        // Recompute map
                        computeNodeToRestrictionMap();
                        updateRestrictedNodePartials = false;
                    }
                    if (partialsMap[nodeNum] != null) {

                    }
                }

            }
        }

        return update;

    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int[] branchUpdateIndices;
    private double[] branchLengths;
    private int branchUpdateCount;

    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;
    //    private final boolean hasRestrictedPartials;
    private static final boolean hasRestrictedPartials = false;

    private final int numRestrictedPartials;
    private final Map<Set<String>, Parameter> partialsRestrictions;
    private Parameter[] partialsMap;
    private double[] partials;
    private boolean updateRestrictedNodePartials;
//    private int[] restrictedIndices;

    protected BufferIndexHelper partialBufferHelper;
    protected BufferIndexHelper scaleBufferHelper;

    protected final int tipCount;
    protected final int internalNodeCount;

    private PartialsRescalingScheme rescalingScheme;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private boolean delayRescalingUntilUnderflow = true;

    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;
//    private int storedRescalingCount;

    private int rescalingMessageCount = 0;

    /**
     * the branch-site model for these sites
     */
    protected final BranchModel branchModel;

    /**
     * A delegate to handle substitution models on branches
     */
    protected final SubstitutionModelDelegate substitutionModelDelegate;

    /**
     * the site model for these sites
     */
    protected final SiteRateModel siteRateModel;

    /**
     * the branch rate model
     */
    protected final BranchRateModel branchRateModel;

    /**
     * the tip partials model
     */
    private final TipStatesModel tipStatesModel;

    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * an array used to transfer tip states
     */
    protected int[] tipStates;

    /**
     * the BEAGLE library instance
     */
    protected Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;

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
    protected final boolean useAmbiguities;

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

            HomogeneousBranchModel homogenousBranchSubstitutionModel = new HomogeneousBranchModel(
                    hky1);

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
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            BranchModel homogeneousBranchModel = new HomogeneousBranchModel(hky1);

            BranchModel epochBranchModel = new EpochBranchModel(treeModel, substitutionModels, epochTimes);

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
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(partitionsList
//            		, sequenceLength
            );
            Alignment alignment = simulator.simulate(false, false);

            BeagleTreeLikelihood nbtl = new BeagleTreeLikelihood(alignment, treeModel, homogeneousBranchModel, siteRateModel, branchRateModel, null, false, PartialsRescalingScheme.DEFAULT, false);

            System.out.println("nBTL(homogeneous) = " + nbtl.getLogLikelihood());

            nbtl = new BeagleTreeLikelihood(alignment, treeModel, epochBranchModel, siteRateModel, branchRateModel, null, false, PartialsRescalingScheme.DEFAULT, false);

            System.out.println("nBTL(epoch) = " + nbtl.getLogLikelihood());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block
    }

    public Double getUpdateTimer() {
        return Double.valueOf(substitutionModelDelegate.updateTime);
    }

    public Double getConvolveTimer() {
        return Double.valueOf(substitutionModelDelegate.convolveTime);
    }

    public void getLogScalingFactors(int nodeIndex, double[] buffer) {
        if (nodeIndex < tipCount) {
            Arrays.fill(buffer, 0.0);
        } else {
//            final int scaleIndex = scaleBufferHelper.getOffsetIndex(nodeIndex - tipCount);
            final int scaleIndex = scaleBufferIndices[nodeIndex - tipCount];
            beagle.getLogScaleFactors(scaleIndex, buffer);
        }
    }

    public double[] getSiteLogLikelihoods() {
        getLogLikelihood();
        double[] siteLogLikelihoods = new double[patternCount];
        beagle.getSiteLogLikelihoods(siteLogLikelihoods);
        return siteLogLikelihoods;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using BEAGLE likelihood calculation library";
    }

    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.AYRES_2012_BEAGLE);
    }

}//END: class