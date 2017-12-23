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

/*
 * ContinuousDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate for continuous traits
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @author Philippe Lemey
 * @version $Id$
 */

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.*;
import dr.evomodel.treedatalikelihood.preorder.ConditionalVarianceAndTransform;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodel.treedatalikelihood.preorder.TipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.TipGradientViaFullConditionalDelegate;
import dr.inference.model.*;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.Matrix;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;
import java.util.logging.Logger;

public class ContinuousDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate,
        ConjugateWishartStatisticsProvider, Citable {

    private final int numTraits;
    private final int dimTrait;
    private final PrecisionType precisionType;
    private final ContinuousRateTransformation rateTransformation;
    private final Tree tree;
    private final BranchRateModel rateModel;
    private final ConjugateRootTraitPrior rootPrior;
    private final boolean forceCompletelyObserved;

    private double branchNormalization;
    private double storedBranchNormalization;

    private boolean allowSingular = false;

    private TreeDataLikelihood callbackLikelihood = null;

    public ContinuousDataLikelihoodDelegate(Tree tree,
                                            DiffusionProcessDelegate diffusionProcessDelegate,
                                            ContinuousTraitPartialsProvider dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            BranchRateModel rateModel,
                                            boolean allowSingular) {
        this(tree, diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel,
                false, allowSingular);
    }

    public ContinuousDataLikelihoodDelegate(Tree tree,
                                            DiffusionProcessDelegate diffusionProcessDelegate,
                                            ContinuousTraitPartialsProvider dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            BranchRateModel rateModel,
                                            boolean forceCompletelyObserved,
                                            boolean allowSingular) {

        super("ContinuousDataLikelihoodDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel.treedatalikelihood");

        logger.info("Using ContinuousDataLikelihood Delegate");

        this.diffusionProcessDelegate = diffusionProcessDelegate;
        addModel(diffusionProcessDelegate);

        this.dataModel = dataModel;
        if (dataModel instanceof Model) {
            addModel((Model) dataModel);
        }

        if (rateModel != null) {
            addModel(rateModel);
        }

        this.numTraits = dataModel.getTraitCount();
        this.dimTrait = dataModel.getTraitDimension();

        this.precisionType = diffusionProcessDelegate.hasDrift() ? // TODO Handle drift in Basic/SCALAR integrator
                PrecisionType.FULL :
                forceCompletelyObserved ?
                        PrecisionType.SCALAR :
                        dataModel.getPrecisionType();
        
        this.rateTransformation = rateTransformation;
        this.tree = tree;
        this.rateModel = rateModel;
        this.rootPrior = rootPrior;

        this.forceCompletelyObserved = forceCompletelyObserved;
        this.allowSingular = allowSingular;

        int nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        int internalNodeCount = nodeCount - tipCount;

        branchUpdateIndices = new int[nodeCount];
        branchLengths = new double[nodeCount];

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

            ContinuousDiffusionIntegrator base;
            if (precisionType == PrecisionType.SCALAR) {

                base = new ContinuousDiffusionIntegrator.Basic(
                        precisionType,
                        numTraits,
                        dimTrait,
                        partialBufferCount,
                        matrixBufferCount
                );

            } else if (precisionType == PrecisionType.FULL) {

                if (diffusionProcessDelegate instanceof DriftDiffusionModelDelegate) {
                    base = new SafeMultivariateWithDriftIntegrator(
                            precisionType,
                            numTraits,
                            dimTrait,
                            partialBufferCount,
                            matrixBufferCount
                    );
                } else {
                    if (allowSingular) {
                        base = new SafeMultivariateIntegrator(
                                precisionType,
                                numTraits,
                                dimTrait,
                                partialBufferCount,
                                matrixBufferCount
                        );
                    } else {
                        base = new MultivariateIntegrator(
                                precisionType,
                                numTraits,
                                dimTrait,
                                partialBufferCount,
                                matrixBufferCount
                        );
                    }
                }

            } else {
                throw new RuntimeException("Not yet implemented");
            }

            cdi = base;
            System.err.println("Base CDI is " + cdi.getClass().getCanonicalName());

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

    TreeDataLikelihood getCallbackLikelihood() { return callbackLikelihood; }

    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    public ContinuousTraitPartialsProvider getDataModel() { return dataModel; }

    public RootProcessDelegate getRootProcessDelegate() {
        return rootProcessDelegate;
    }

    public ConjugateRootTraitPrior getRootPrior() { return rootPrior; }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();
        sb.append("Tree:\n");
        sb.append(callbackLikelihood.getId()).append("\t");

        sb.append(cdi.getReport());

        final Tree tree = callbackLikelihood.getTree();
        sb.append(tree.toString());
        sb.append("\n\n");

        final double normalization = rateTransformation.getNormalization();
        final double priorSampleSize = rootProcessDelegate.getPseudoObservations();

        double[][] treeStructure = MultivariateTraitDebugUtilities.getTreeVariance(tree, 1.0, Double.POSITIVE_INFINITY);
        sb.append("Tree structure:\n");
        sb.append(new Matrix(treeStructure));
        sb.append("\n\n");

        double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree, normalization, priorSampleSize);
        double[][] traitPrecision = getDiffusionModel().getPrecisionmatrix();
        Matrix traitVariance = new Matrix(traitPrecision).inverse();

        double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());

        Matrix treeV = new Matrix(treeVariance);
        Matrix treeP = treeV.inverse();

        sb.append("Tree variance:\n");
        sb.append(treeV);
        sb.append("Tree precision:\n");
        sb.append(treeP);
        sb.append("\n\n");

        sb.append("Trait variance:\n");
        sb.append(traitVariance);
        sb.append("\n\n");

        sb.append("Joint variance:\n");
        sb.append(new Matrix(jointVariance));
        sb.append("\n\n");

        double[][] treeDrift = MultivariateTraitDebugUtilities.getTreeDrift(tree, diffusionProcessDelegate);
        if (diffusionProcessDelegate.hasDrift()) {
            sb.append("Tree drift:\n");
            sb.append(new Matrix(treeDrift));
            sb.append("\n\n");
        }
        double[] drift = KroneckerOperation.vectorize(treeDrift);

        final int datumLength = tipCount * dimTrait;

        sb.append("Tree dim : ").append(treeVariance.length).append("\n");
        sb.append("dimTrait : ").append(dimTrait).append("\n");
        sb.append("numTraits: ").append(numTraits).append("\n");
        sb.append("jVar dim : ").append(jointVariance.length).append("\n");
        sb.append("datum dim: ").append(datumLength);
        sb.append("\n\n");

        double[] data = dataModel.getParameter().getParameterValues();

        if (dataModel instanceof ContinuousTraitDataModel) {
            for (int tip = 0; tip < tipCount; ++tip) {
                double[] tipData = ((ContinuousTraitDataModel) dataModel).getTipObservation(tip, precisionType);
                System.arraycopy(tipData, 0, data, tip * numTraits * dimTrait, numTraits * dimTrait);
            }
        }

        sb.append("data: ").append(new dr.math.matrixAlgebra.Vector(data));
        sb.append("\n\n");

        double[][] graphStructure = MultivariateTraitDebugUtilities.getGraphVariance(tree, 1.0,
                priorSampleSize);
        double[][] jointGraphVariance = KroneckerOperation.product(graphStructure, traitVariance.toComponents());

        sb.append("graph structure:\n");
        sb.append(new Matrix(graphStructure));
        sb.append("\n\n");

        double[] priorMean = rootPrior.getMean();
        sb.append("prior mean: ").append(new dr.math.matrixAlgebra.Vector(priorMean));
        sb.append("\n\n");

        for (int index = 0; index < drift.length / dimTrait; ++index) {
            for (int dim = 0; dim < dimTrait; ++dim) {
                drift[index * dimTrait + dim] += priorMean[dim];
            }
        }

        for (int trait = 0; trait < numTraits; ++trait) {
            sb.append("Trait #").append(trait).append("\n");

            double[] rawDatum = new double[datumLength];

            List<Integer> missing = new ArrayList<Integer>();
            int index = 0;
            for (int tip = 0; tip < tipCount; ++tip) {
                for (int dim = 0; dim < dimTrait; ++dim) {
                    double d = data[tip * dimTrait * numTraits + trait * dimTrait + dim];
                    rawDatum[index] = d;
                    if (Double.isNaN(d)) {
                        missing.add(index);
                    }
                    ++index;
                }
            }

            double[][] varianceDatum = jointVariance;
            double[] datum = rawDatum;

            double[] driftDatum = drift;
            
            int[] notMissingIndices;
            if (missing.size() > 0) {
                notMissingIndices = new int[datumLength - missing.size()];
                int offsetNotMissing = 0;
                for (int i = 0; i < datumLength; ++i) {
                    if (!missing.contains(i)) {
                        notMissingIndices[offsetNotMissing] = i;
                        ++offsetNotMissing;
                    }
                }

                datum = Matrix.gatherEntries(rawDatum, notMissingIndices);
                varianceDatum = Matrix.gatherRowsAndColumns(jointVariance, notMissingIndices, notMissingIndices);

                driftDatum = Matrix.gatherEntries(drift, notMissingIndices);
            }

            sb.append("datum : ").append(new dr.math.matrixAlgebra.Vector(datum)).append("\n");

            sb.append("drift : ").append(new dr.math.matrixAlgebra.Vector(driftDatum)).append("\n");

            sb.append("variance:\n");
            sb.append(new Matrix(varianceDatum));

            MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(driftDatum, new Matrix(varianceDatum).inverse().toComponents());
            double logDensity = mvn.logPdf(datum);
            sb.append("\n\n");
            sb.append("logDatumLikelihood: ").append(logDensity).append("\n\n");

            // Compute joint for internal nodes
            int[] cNotMissingJoint = new int[dimTrait * tipCount];
            int[] cMissingJoint = new int[dimTrait * (tipCount - 1)];

            // External nodes
            for (int tipTrait = 0; tipTrait < dimTrait * tipCount; ++tipTrait) {
                cNotMissingJoint[tipTrait] = tipTrait;
            }

            // Internal nodes
            for (int tipTrait = dimTrait * tipCount; tipTrait < dimTrait * (2 * tipCount - 1); ++tipTrait) {
                cMissingJoint[tipTrait - dimTrait * tipCount] = tipTrait;
            }
            
            double[] rawDatumJoint = new double[dimTrait * (2 * tipCount - 1)];
            System.arraycopy(rawDatum, 0, rawDatumJoint, 0, rawDatum.length);

            double[][] driftJointMatrix = MultivariateTraitDebugUtilities.getGraphDrift(tree, diffusionProcessDelegate);
            double[] driftJoint = KroneckerOperation.vectorize(driftJointMatrix);
            
            for (int idx = 0; idx < driftJoint.length / dimTrait; ++idx) {
                for (int dim = 0; dim < dimTrait; ++dim) {
                    driftJoint[idx * dimTrait + dim] += priorMean[dim];
                }
            }

            ConditionalVarianceAndTransform cVarianceJoint = new ConditionalVarianceAndTransform(
                    new Matrix(jointGraphVariance), cMissingJoint, cNotMissingJoint);

            double[] cMeanJoint = cVarianceJoint.getConditionalMean(rawDatumJoint, 0, driftJoint, 0);

            sb.append("cDriftJoint: ").append(new dr.math.matrixAlgebra.Vector(driftJoint)).append("\n\n");

            sb.append("cMeanInternalJoint: ").append(new dr.math.matrixAlgebra.Vector(cMeanJoint)).append("\n\n");

            // Compute full conditional distributions
            sb.append("Full conditional distributions:\n");


            int[] cMissing = new int[dimTrait];
            int[] cNotMissing = new int[tipCount * dimTrait - dimTrait];

            for (int tip = 0; tip < tipCount; ++tip) {

                for (int cTrait = 0; cTrait < dimTrait; ++cTrait) {
                    cMissing[cTrait] = tip * dimTrait + cTrait;
                }

                for (int m = 0; m < tip * dimTrait; ++m) {
                    cNotMissing[m] = m;
                }

                for (int m = (tip + 1) * dimTrait; m < tipCount * dimTrait; ++m) {
                    cNotMissing[m - dimTrait] = m;
                }

                ConditionalVarianceAndTransform cVariance = new ConditionalVarianceAndTransform(
                        new Matrix(jointVariance), cMissing, cNotMissing);

                double[] cMean = cVariance.getConditionalMean(rawDatum, 0, drift, 0);

                sb.append("cMean #").append(tip).append(" ").append(new dr.math.matrixAlgebra.Vector(cMean))
                    .append("\n");
            }
        }

        return sb.toString();
    }

    @Override
    public final int getTraitCount() {
        return numTraits;
    }

    @Override
    public final int getTraitDim() {
        return dimTrait;
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        return rateTransformation.getRateRescalingScheme();
    }

    public final ContinuousDiffusionIntegrator getIntegrator() {
        return cdi;
    }

    final ContinuousRateTransformation getRateTransformation() {
        return rateTransformation;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
        this.callbackLikelihood = treeDataLikelihood;
    }

    @Override
    public int vectorizeNodeOperations(final List<ProcessOnTreeDelegate.NodeOperation> nodeOperations,
                                       final int[] operations) {

        int k = 0;
        for (NodeOperation op : nodeOperations) {

            operations[k    ] = getActiveNodeIndex(op.getNodeNumber());
            operations[k + 1] = getActiveNodeIndex(op.getLeftChild());    // source node 1
            operations[k + 2] = getActiveMatrixIndex(op.getLeftChild());  // source matrix 1
            operations[k + 3] = getActiveNodeIndex(op.getRightChild());   // source node 2
            operations[k + 4] = getActiveMatrixIndex(op.getRightChild()); // source matrix 2

            k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
        }

        return nodeOperations.size();
    }

//    public static String getStringOfVectorizedOperations(final int[] operations, final int count) {
//        StringBuilder sb = new StringBuilder();
//        int k = 0;
//        for (int i = 0; i < count; ++i) {
//            sb.append(operations[k    ]).append(" ");
//            sb.append(operations[k + 1]).append(" ");
//            sb.append(operations[k + 2]).append(" ");
//            sb.append(operations[k + 3]).append(" ");
//            sb.append(operations[k + 4]).append("\n");
//
//            k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
//        }
//
//        return sb.toString();
//    }

    public DiffusionProcessDelegate getDiffusionProcessDelegate() { return diffusionProcessDelegate; }

    public MultivariateDiffusionModel getDiffusionModel() {
        return diffusionProcessDelegate.getDiffusionModel(0);
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

        final double[] tipPartial = dataModel.getTipPartial(tipIndex, forceCompletelyObserved);
        setTipDataDirectly(tipIndex, tipPartial);
    }

    void setTipDataDirectly(int tipIndex, double[] tipPartial) {
        cdi.setPostOrderPartial(partialBufferHelper.getOffsetIndex(tipIndex),
                tipPartial);
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
                                      int rootNodeNumber) throws LikelihoodException {

        branchNormalization = rateTransformation.getNormalization();

        int branchUpdateCount = 0;
        for (BranchOperation op : branchOperations) {
            branchUpdateIndices[branchUpdateCount] = op.getBranchNumber();
            branchLengths[branchUpdateCount] = op.getBranchLength() * branchNormalization;
            branchUpdateCount ++;
        }

        if (!updateTipData.isEmpty()) {
            if (updateTipData.getFirst() == -1) { // Update all tips
                setAllTipData(flip);
            } else {
                while(!updateTipData.isEmpty()) {
                    int tipIndex = updateTipData.removeFirst();
                    setTipData(tipIndex, flip);
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

        int operationCount = vectorizeNodeOperations(nodeOperations, operations);

        int[] degreesOfFreedom = null;
        double[] outerProducts = null;

        if (computeWishartStatistics) {
            degreesOfFreedom = new int[numTraits];
            outerProducts = new double[dimTrait * dimTrait  * numTraits];
            cdi.setWishartStatistics(degreesOfFreedom, outerProducts);
        }

        cdi.updatePostOrderPartials(operations, operationCount, computeWishartStatistics);

        double[] logLikelihoods = new double[numTraits];

        rootProcessDelegate.calculateRootLogLikelihood(cdi, partialBufferHelper.getOffsetIndex(rootNodeNumber),
                logLikelihoods, computeWishartStatistics);

        if (computeWishartStatistics) {
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

    public final int getActiveNodeIndex(final int index) {
        return partialBufferHelper.getOffsetIndex(index);
    }

    public final int getActiveMatrixIndex(final int index) {
        return diffusionProcessDelegate.getMatrixIndex(index);
    }

    public void getPostOrderPartial(final int nodeNumber, double[] vector) {
        cdi.getPostOrderPartial(getActiveNodeIndex(nodeNumber), vector);
    }

//    public void getPostOrderPartial(final int nodeNumber, double[] vector, double[] matrix, double[] displacement) {
//        cdi.getPostOrderPartial(getActiveNodeIndex(nodeNumber), vector, matrix, displacement);
//    }
//
//    public void getPreOrderPartial(final int nodeNumber, double[] vector) {
//        cdi.getPreOrderPartial(getActiveNodeIndex(nodeNumber), vector);
//    }

    @Override
    public void makeDirty() {
        updateDiffusionModel = true;
        fireModelChanged(); // Signal simulation processes
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == diffusionProcessDelegate) {
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

        // turn on double buffering flipping (may have been turned off to enable a rescale)
        flip = true;

        storedBranchNormalization = branchNormalization;
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restoreState() {
        partialBufferHelper.restoreState();

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

    private final int tipCount;

    private final int[] branchUpdateIndices;
    private final double[] branchLengths;

    private final int[] operations;

    private boolean flip = true;
    private final BufferIndexHelper partialBufferHelper;

    private final DiffusionProcessDelegate diffusionProcessDelegate;

    private final RootProcessDelegate rootProcessDelegate;

    private final ContinuousTraitPartialsProvider dataModel;

    private final ContinuousDiffusionIntegrator cdi;

    private boolean updateDiffusionModel;

    private final Deque<Integer> updateTipData = new ArrayDeque<Integer>();

    private WishartSufficientStatistics wishartStatistics = null;

    private boolean computeWishartStatistics = false;

    @Override
    public WishartSufficientStatistics getWishartStatistics() {
        return wishartStatistics;
    }

    void setComputeWishartStatistics(boolean computeWishartStatistics) {
        this.computeWishartStatistics = computeWishartStatistics;
    }

    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return getDiffusionModel().getPrecisionParameter();
    }

    public void addFullConditionalGradientTrait(String traitName) {

        ProcessSimulationDelegate gradientDelegate = new TipGradientViaFullConditionalDelegate(traitName,
                (MutableTreeModel) getCallbackLikelihood().getTree(),
                getDiffusionModel(),
                (ContinuousTraitDataModel) getDataModel(), getRootPrior(),
                getRateTransformation(), this);

        TreeTraitProvider traitProvider = new ProcessSimulation(getCallbackLikelihood(), gradientDelegate);

        getCallbackLikelihood().addTraits(traitProvider.getTreeTraits());
    }

    public void addFullConditionalDensityTrait(String traitName) {

        ProcessSimulationDelegate gradientDelegate = new TipFullConditionalDistributionDelegate(traitName,
                getCallbackLikelihood().getTree(),
                getDiffusionModel(),
                getDataModel(), getRootPrior(),
                getRateTransformation(), this);

        TreeTraitProvider traitProvider = new ProcessSimulation(getCallbackLikelihood(), gradientDelegate);

        getCallbackLikelihood().addTraits(traitProvider.getTreeTraits());
    }

    static ContinuousDataLikelihoodDelegate createObservedDataOnly(ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        return new ContinuousDataLikelihoodDelegate(likelihoodDelegate.tree,
                likelihoodDelegate.diffusionProcessDelegate,
                likelihoodDelegate.dataModel,
                likelihoodDelegate.rootPrior,
                likelihoodDelegate.rateTransformation,
                likelihoodDelegate.rateModel,
                true,
                likelihoodDelegate.allowSingular);
    }

    static ContinuousDataLikelihoodDelegate createWithMissingData(ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        if (!(likelihoodDelegate.dataModel instanceof ContinuousTraitDataModel)) {
            throw new IllegalArgumentException("Not yet implemented");
        }

        List<Integer> originalMissingIndices = ((ContinuousTraitDataModel) likelihoodDelegate.dataModel).getOriginalMissingIndices();

        if (originalMissingIndices.size() == 0) {
            throw new IllegalArgumentException("ContinuousDataLikelihoodDelegate has no missing traits");
        }

        ContinuousTraitPartialsProvider newDataModel = new ContinuousTraitDataModel(((ContinuousTraitDataModel) likelihoodDelegate.dataModel).getName(),
                likelihoodDelegate.dataModel.getParameter(),
                ((ContinuousTraitDataModel) likelihoodDelegate.dataModel).getOriginalMissingIndices(),
                true,
                likelihoodDelegate.getTraitDim(), PrecisionType.FULL);

        return new ContinuousDataLikelihoodDelegate(likelihoodDelegate.tree,
                likelihoodDelegate.diffusionProcessDelegate,
                newDataModel,
                likelihoodDelegate.rootPrior,
                likelihoodDelegate.rateTransformation,
                likelihoodDelegate.rateModel,
                false,
                likelihoodDelegate.allowSingular);
    }
}
