/*
 * ContinuousDataLikelihoodDelegate.java
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

package dr.evomodel.treedatalikelihood.continuous;

/**
 * ContinuousDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate for continuous traits
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @author Philippe Lemey
 * @version $Id$
 */

import com.sun.org.apache.regexp.internal.RE;
import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;
import java.util.logging.Logger;

public class ContinuousDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate, ConjugateWishartStatisticsProvider, Citable {

    private final int numTraits;
    private final int dimTrait;
    private final PrecisionType precisionType;
    private final ContinuousRateTransformation rateTransformation;

    private double branchNormalization;
    private double storedBranchNormalization;

    private TreeDataLikelihood callbackLikelihood = null;

    public ContinuousDataLikelihoodDelegate(MultivariateTraitTree tree,
                                            MultivariateDiffusionModel diffusionModel,
                                            ContinuousTraitDataModel dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            BranchRateModel rateModel) {

        super("ContinousDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel.treedatalikelihood");

        logger.info("Using ContinuousDataLikelihood Delegate");

        this.diffusionModel = diffusionModel;
        addModel(diffusionModel);

        this.dataModel = dataModel;
        addModel(dataModel);

        if (rateModel != null) {
            addModel(rateModel);
        }

        this.numTraits = dataModel.getTraitCount();
        this.dimTrait = dataModel.getTraitDimension();
        this.precisionType = dataModel.getPrecisionType();
        this.rateTransformation = rateTransformation;

        nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        branchUpdateIndices = new int[nodeCount];
        branchLengths = new double[nodeCount];

        diffusionProcessDelegate = new HomogenousDiffusionModelDelegate(tree, diffusionModel);

        // one or two partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(nodeCount,
                dataModel.bufferTips() ? 0 : tipCount);

        int partialBufferCount = partialBufferHelper.getBufferCount();
        int matrixBufferCount = diffusionProcessDelegate.getEigenBufferCount();

        rootProcessDelegate = new RootProcessDelegate.FullyConjugate(rootPrior, precisionType, numTraits,
                partialBufferCount, matrixBufferCount);

        partialBufferCount += rootProcessDelegate.getExtraPartialBufferCount();
        matrixBufferCount += rootProcessDelegate.getExtraMatrixBufferCount();

        operations = new int[(internalNodeCount + rootProcessDelegate.getExtraPartialBufferCount())
                * ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE];

        try {

            boolean USE_OLD = false;

            if (precisionType == PrecisionType.SCALAR || USE_OLD) {

                cdi = new ContinuousDiffusionIntegrator.Basic(
                        precisionType,
                        numTraits,
                        dimTrait,
                        partialBufferCount,
                        matrixBufferCount
                );

            } else if (precisionType == PrecisionType.FULL) {

                cdi = new ContinuousDiffusionIntegrator.Multivariate(
                        precisionType,
                        numTraits,
                        dimTrait,
                        partialBufferCount,
                        matrixBufferCount
                );

            } else {
                throw new RuntimeException("Not yet implemented");
            }

            // TODO Make separate library
//            cdi = CDIFactory.loadCDIInstance();
//
//            InstanceDetails instanceDetails = cdi.getDetails();
//            ResourceDetails resourceDetails = null;
//
//            if (instanceDetails != null) {
//                resourceDetails = CDIFactory.getResourceDetails(instanceDetails.getResourceNumber());
//                if (resourceDetails != null) {
//                    StringBuilder sb = new StringBuilder("  Using CDI resource ");
//                    sb.append(resourceDetails.getNumber()).append(": ");
//                    sb.append(resourceDetails.getName()).append("\n");
//                    if (resourceDetails.getDescription() != null) {
//                        String[] description = resourceDetails.getDescription().split("\\|");
//                        for (String desc : description) {
//                            if (desc.trim().length() > 0) {
//                                sb.append("    ").append(desc.trim()).append("\n");
//                            }
//                        }
//                    }
//                    sb.append("    with instance flags: ").append(instanceDetails.toString());
//                    logger.info(sb.toString());
//                } else {
//                    logger.info("  Error retrieving CDI resource for instance: " + instanceDetails.toString());
//                }
//            } else {
//                logger.info("  No external CDI resources available, or resource list/requirements not met, using Java implementation");
//            }

            // Check tip data
            for (int i = 0; i < tipCount; i++) {
                final NodeRef node = tree.getExternalNode(i);
                final int index = node.getNumber();
                
                assert (i == index);

                if (!checkDataAlignment(node, tree)) {
                    throw new TaxonList.MissingTaxonException("Missing taxon");
                }
            }            

            setAllTipData(dataModel.bufferTips());

            rootProcessDelegate.setRootPartial(cdi);

            updateDiffusionModel = true;
            
        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
    }

    @Override
    public final int getTraitCount() {
        return numTraits;
    }

    @Override
    public final int getTraitDim() {
        return dimTrait;
    }

    public final ContinuousDiffusionIntegrator getIntegrator() {
        return cdi;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
        this.callbackLikelihood = treeDataLikelihood;
    }

    public MultivariateDiffusionModel getDiffusionModel() {
        return diffusionModel;
    }

    private void setAllTipData(boolean flip) {
        for (int index = 0; index < tipCount; index++) {
             setTipData(index, flip);
        }
        updateTipData.clear();
    }
    
    private void setTipData(int tipIndex, boolean flip) {
        if (flip) {
            partialBufferHelper.flipOffset(tipIndex);
        }

        cdi.setPartial(partialBufferHelper.getOffsetIndex(tipIndex),
                dataModel.getTipPartial(tipIndex));
    }
    
    private boolean checkDataAlignment(NodeRef node, Tree tree) {
        int index = node.getNumber();
        String name1 = dataModel.getParameter().getParameter(index).getParameterName();
        Taxon taxon = tree.getNodeTaxon(node);
        return name1.contains(taxon.getId());
    }

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        return TreeTraversal.TraversalType.POST_ORDER;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations,
                                      int rootNodeNumber) throws LikelihoodUnderflowException {

        // TODO REMOVE NEXT LINE
        computeWishartStatistics = true;

        branchNormalization = rateTransformation.getNormalization();  // TODO Cache branchNormalization

        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchLength() * branchNormalization;
            branchUpdateCount ++;
        }

        if (!updateTipData.isEmpty()) {
            if (updateTipData.getFirst() == -1) { // Update all tips
                setAllTipData(flip);
//                System.err.println("SET ALL TIPS");
            } else {
                while(!updateTipData.isEmpty()) {
                    int tipIndex = updateTipData.removeFirst();
                    setTipData(tipIndex, flip);
//                    System.err.println("SET TIP " + tipIndex);
                }
            }
        }

        if (updateDiffusionModel) {
            diffusionProcessDelegate.setDiffusionModels(cdi, flip);
        }

        if (branchUpdateCount > 0) {
            diffusionProcessDelegate.updateDiffusionMatrices(
                    cdi,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount,
                    flip);
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

            operations[k + 0] = partialBufferHelper.getOffsetIndex(nodeNum);
            operations[k + 1] = partialBufferHelper.getOffsetIndex(op.getLeftChild()); // source node 1
            operations[k + 2] = diffusionProcessDelegate.getMatrixIndex(op.getLeftChild()); // source matrix 1
            operations[k + 3] = partialBufferHelper.getOffsetIndex(op.getRightChild()); // source node 2
            operations[k + 4] = diffusionProcessDelegate.getMatrixIndex(op.getRightChild()); // source matrix 2

            k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
        }

        int[] degreesOfFreedom = null;
        double[] outerProducts = null;

        if (computeWishartStatistics) {
            // TODO Abstract this ugliness away
            degreesOfFreedom = new int[numTraits];
            outerProducts = new double[dimTrait * dimTrait  * numTraits];
            cdi.setWishartStatistics(degreesOfFreedom, outerProducts);
        }

        cdi.updatePartials(operations, operationCount, computeWishartStatistics);

        double[] logLikelihoods = new double[numTraits];

        rootProcessDelegate.calculateRootLogLikelihood(cdi, partialBufferHelper.getOffsetIndex(rootNodeNumber),
                logLikelihoods, computeWishartStatistics);

        if (computeWishartStatistics) {
            // TODO Abstract this ugliness away
            cdi.getWishartStatistics(degreesOfFreedom, outerProducts);
            wishartStatistics = new WishartSufficientStatistics(
                    degreesOfFreedom,
                    outerProducts
            );

        } else {
            wishartStatistics = null;
        }

        double logL = 0.0;
        for (double d : logLikelihoods) {
            logL += d;
        }

        updateDiffusionModel = false;

        return logL;
    }

    public void getPartial(final int nodeNumber, double[] vector) {
        cdi.getPartial(partialBufferHelper.getOffsetIndex(nodeNumber), vector);
    }

    @Override
    public void makeDirty() {
        updateDiffusionModel = true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == diffusionModel) {
            updateDiffusionModel = true;
            // Tell TreeDataLikelihood to update all nodes
            fireModelChanged();
        } else if (model == dataModel) {
            if (object == dataModel) {
                if (index == -1) { // all taxa updated
                    updateTipData.addFirst(index);
                    fireModelChanged();
                } else { // only one taxon updated
                    updateTipData.addLast(index);
                    fireModelChanged(this, index);
                }
            }

        } else if (model instanceof BranchRateModel) {
            fireModelChanged();
        } else {
            throw new RuntimeException("Unknown model component");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void storeState() {
        partialBufferHelper.storeState();
        diffusionProcessDelegate.storeState();

        // turn on double buffering flipping (may have been turned off to enable a rescale)
        flip = true;

        storedBranchNormalization = branchNormalization;
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restoreState() {
//        updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore

        partialBufferHelper.restoreState();
        diffusionProcessDelegate.restoreState();

        branchNormalization = storedBranchNormalization;
    }

    @Override
    protected void acceptState() {
    }

    // **************************************************************
    // INSTANCE CITABLE
    // **************************************************************

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "TODO";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2010_PHYLOGEOGRAPHY);
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final int nodeCount;
    private final int tipCount;
    private final int internalNodeCount;

    private final int[] branchUpdateIndices;
    private final double[] branchLengths;

    private final int[] operations;

    private boolean flip = true;
    private final BufferIndexHelper partialBufferHelper;

    private final DiffusionProcessDelegate diffusionProcessDelegate;

    private final MultivariateDiffusionModel diffusionModel;

    private final RootProcessDelegate rootProcessDelegate;

    private final ContinuousTraitDataModel dataModel;

    private final ContinuousDiffusionIntegrator cdi;

    private boolean updateDiffusionModel;

    private final Deque<Integer> updateTipData = new ArrayDeque<Integer>();

    private WishartSufficientStatistics wishartStatistics = null;

    private boolean computeWishartStatistics = false;

    @Override
    public WishartSufficientStatistics getWishartStatistics() {
        assert (callbackLikelihood != null);
        callbackLikelihood.makeDirty();
        computeWishartStatistics = true;
        callbackLikelihood.getLogLikelihood();
        computeWishartStatistics = false;
        return wishartStatistics;
    }

    @Override
    public MatrixParameterInterface getPrecisionParamter() {
        return diffusionModel.getPrecisionParameter();
    }

//    private boolean updateSiteModel;
}
