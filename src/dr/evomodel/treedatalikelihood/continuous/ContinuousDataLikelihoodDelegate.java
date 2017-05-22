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
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.WrappedVector;
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
    private final MultivariateTraitTree tree;
    private final BranchRateModel rateModel;
    private final ConjugateRootTraitPrior rootPrior;
    private final boolean forceCompletelyObserved;

    private double branchNormalization;
    private double storedBranchNormalization;


    private TreeDataLikelihood callbackLikelihood = null;

    public ContinuousDataLikelihoodDelegate(MultivariateTraitTree tree,
                                            MultivariateDiffusionModel diffusionModel,
                                            ContinuousTraitDataModel dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            BranchRateModel rateModel) {
        this(tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, false);
    }

    public ContinuousDataLikelihoodDelegate(MultivariateTraitTree tree,
                                            MultivariateDiffusionModel diffusionModel,
                                            ContinuousTraitDataModel dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            BranchRateModel rateModel,
                                            boolean forceCompletelyObserved) {

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
        this.precisionType = forceCompletelyObserved ? PrecisionType.SCALAR : dataModel.getPrecisionType();
        this.rateTransformation = rateTransformation;
        this.tree = tree;
        this.rateModel = rateModel;
        this.rootPrior = rootPrior;
        this.forceCompletelyObserved = forceCompletelyObserved;

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

            ContinuousDiffusionIntegrator base = null;

            if (precisionType == PrecisionType.SCALAR || USE_OLD) {

                base = new ContinuousDiffusionIntegrator.Basic(
                        precisionType,
                        numTraits,
                        dimTrait,
                        partialBufferCount,
                        matrixBufferCount
                );

            } else if (precisionType == PrecisionType.FULL) {

                base = new ContinuousDiffusionIntegrator.Multivariate(
                        precisionType,
                        numTraits,
                        dimTrait,
                        partialBufferCount,
                        matrixBufferCount
                );

            } else {
                throw new RuntimeException("Not yet implemented");
            }

//            cdi = new ContinuousDiffusionIntegrator.OuterProductProvider(base);
            cdi = base;

            System.err.println("Base CDI is " + cdi.getClass().getCanonicalName());
//            System.exit(-1);

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

    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    private double[] getTipObservations() {
        final double[] data = new double[numTraits * dimTrait * tipCount];

        for (int tip = 0; tip < tipCount; ++tip) {
            double[] tipData = dataModel.getTipObservation(tip, precisionType);
            System.arraycopy(tipData, 0, data, tip * numTraits * dimTrait, numTraits * dimTrait);
        }

        return data;
    }

    public RootProcessDelegate getRootProcessDelegate() {
        return rootProcessDelegate;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();
        sb.append("Tree:\n");
        sb.append(callbackLikelihood.getId()).append("\t");

        sb.append(cdi.getReport());

//        final Tree tree = callbackLikelihood.getTree();
//        sb.append(tree.toString());
//        sb.append("\n\n");
//
//        final double normalization = rateTransformation.getNormalization();
//        final double priorSampleSize = rootProcessDelegate.getPseudoObservations();
//
//        double[][] treeStructure = MultivariateTraitDebugUtilities.getTreeVariance(tree, 1.0, Double.POSITIVE_INFINITY);
//        sb.append("Tree structure:\n");
//        sb.append(new Matrix(treeStructure));
//        sb.append("\n\n");
//
//        double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree, normalization, priorSampleSize);
//        double[][] traitPrecision = getDiffusionModel().getPrecisionmatrix();
//        Matrix traitVariance = new Matrix(traitPrecision).inverse();
//
//        double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());
//
//        Matrix treeV = new Matrix(treeVariance);
//        Matrix treeP = treeV.inverse();
//
//        sb.append("Tree variance:\n");
//        sb.append(treeV);
//        sb.append("Tree precision:\n");
//        sb.append(treeP);
////        sb.append(matrixMin(treeVariance)).append("\t").append(matrixMax(treeVariance)).append("\t").append(matrixSum(treeVariance));
//        sb.append("\n\n");
//        sb.append("Trait variance:\n");
//        sb.append(traitVariance);
//        sb.append("\n\n");
//        sb.append("Joint variance:\n");
//        sb.append(new Matrix(jointVariance));
//        sb.append("\n\n");
//
//        final int datumLength = tipCount * dimTrait;
//
//        sb.append("Tree dim : " + treeVariance.length + "\n");
//        sb.append("dimTrait : " + dimTrait + "\n");
//        sb.append("numTraits: " + numTraits + "\n");
//        sb.append("Jvar dim : " + jointVariance.length + "\n");
//        sb.append("datum dim: " + datumLength);
//        sb.append("\n\n");
//
//        double[] data = getTipObservations();
//        sb.append("data: " + new dr.math.matrixAlgebra.Vector(data));
//        sb.append("\n\n");
//
//        double logLikelihood = 0;
//
//        Matrix totalNop = new Matrix(dimTrait, dimTrait);
//        Matrix totalOp = new Matrix(dimTrait, dimTrait);
//
//        for (int trait = 0; trait < numTraits; ++trait) {
//            sb.append("Trait #" + trait + "\n");
//
//            double[] rawDatum = new double[datumLength];
//            double[][] opDatum = new double[tipCount][dimTrait];
//
//
//            List<Integer> missing = new ArrayList<Integer>();
//            int index = 0;
//            for (int tip = 0; tip < tipCount; ++tip) {
//                for (int dim = 0; dim < dimTrait; ++dim) {
//                    double d = data[tip * dimTrait * numTraits + trait * dimTrait + dim];
//                    rawDatum[index] = d;
//                    if (Double.isNaN(d)) {
//                        missing.add(index);
//                        d = 0.0;
//                    }
//                    opDatum[tip][dim] = d;
//                    ++index;
//                }
//            }
//
//            double[][] varianceDatum = jointVariance;
//            double[] datum = rawDatum;
//
//            int[] missingIndices = null;
//            int[] notMissingIndices = null;
//
//            if (missing.size() > 0) {
//                missingIndices = new int[missing.size()];
//                notMissingIndices = new int[datumLength - missing.size()];
//                int offsetMissing = 0;
//                int offsetNotMissing = 0;
//                for (int i = 0; i < datumLength; ++i) {
//                    if (!missing.contains(i)) {
//                        notMissingIndices[offsetNotMissing] = i;
//                        ++offsetNotMissing;
//                    } else {
//                        missingIndices[offsetMissing] = i;
//                        ++offsetMissing;
//                    }
//                }
//
//                datum = Matrix.gatherEntries(rawDatum, notMissingIndices);
//                varianceDatum = Matrix.gatherRowsAndColumns(jointVariance, notMissingIndices, notMissingIndices);
//            }
//
//            sb.append("datum : " + new dr.math.matrixAlgebra.Vector(datum) + "\n");
//            sb.append("variance:\n");
//            sb.append(new Matrix(varianceDatum));
//
//            MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new double[datum.length], new Matrix(varianceDatum).inverse().toComponents());
//            double logDensity = mvn.logPdf(datum);
//            sb.append("\n\n");
//            sb.append("logDatumLikelihood: " + logDensity + "\n\n");
//            logLikelihood += logDensity;
//
//            if (DEBUG_MISSING_DISTRIBUTION && missing.size() > 0) {
//                sb.append("\nConditional distribution of missing values at");
//                for (int m : missing) {
//                    sb.append(" " + m);
//                }
//                sb.append("\n");
////                for (int n : notMissingIndices) {
////                    sb.append(" " + n);
////                }
////                sb.append("\n");
//
//
//                ProcessSimulationDelegate.ConditionalOnPartiallyMissingTipsDelegate.ConditionalVarianceAndTranform transform =
//                        new ProcessSimulationDelegate.ConditionalOnPartiallyMissingTipsDelegate.ConditionalVarianceAndTranform(
//                        new Matrix(jointVariance), missingIndices, notMissingIndices
//                );
//
//                double[] mean = transform.getConditionalMean(rawDatum, 0, new double[rawDatum.length], 0);
//                Matrix variance = transform.getVariance();
//
//                sb.append("obs: " + new WrappedVector.Raw(rawDatum, 0, rawDatum.length));
//                sb.append("cMean: " + new dr.math.matrixAlgebra.Vector(mean) + "\n");
//                sb.append("cVar :\n" + variance + "\n");
////                System.err.println(sb.toString());
////                System.exit(-1);
//            }
//
//            if (DEBUG_OUTER_PRODUCTS) {
//
//                Matrix y = new Matrix(opDatum);
//
////            System.err.println("y = \n" + y);
//                sb.append("Y:\n" + y);
//                sb.append("Tree V:\n" + treeV);
//
//                Matrix op = null;
//
//                try {
//                    op = y.transpose().product(treeP).product(y);
//                    totalOp.accumulate(op);
//                } catch (IllegalDimension illegalDimension) {
//                    illegalDimension.printStackTrace();
//                }
//
//                sb.append("Outer-products:\n");
//                sb.append(op);
//                sb.append("\n\n");
//
//                sb.append("check for missing taxa ...");
//
//                missing.clear();
//                for (int tip = 0; tip < tipCount; ++tip) {
//                    if (allZero(opDatum[tip])) {
//                        missing.add(tip);
//                    }
//                }
//
//                index = 0;
//                int[] notMissing = new int[opDatum.length - missing.size()];
//                double[][] nopDatum = new double[opDatum.length - missing.size()][];
//                for (int tip = 0; tip < tipCount; ++tip) {
//                    if (!missing.contains(tip)) {
//                        nopDatum[index] = opDatum[tip];
//                        notMissing[index] = tip;
//                        ++index;
//                    }
//                }
//
//                Matrix nonMissingTreeVariance = treeV.extractRowsAndColumns(notMissing, notMissing);
//                Matrix notMissingTreePrecision = nonMissingTreeVariance.inverse();
//                Matrix notMissingY = new Matrix(nopDatum);
//
//                sb.append("NP Y:\n" + notMissingY);
//                sb.append("NP Tree V:\n" + nonMissingTreeVariance);
//                sb.append("NP Tree P:\n" + notMissingTreePrecision);
//
//                Matrix nop = null;
//                try {
//                    nop = notMissingY.transpose().product(notMissingTreePrecision).product(notMissingY);
//                    totalNop.accumulate(nop);
//                } catch (IllegalDimension illegalDimension) {
//                    illegalDimension.printStackTrace();
//                }
//
//                sb.append("NP Outer-products:\n");
//                sb.append(nop);
//            }
//
//            sb.append("\n\n");
//
//
//        }
//
//        sb.append("TOTAL DEBUG logLikelihood = " + logLikelihood + "\n");
//
//        if (DEBUG_OUTER_PRODUCTS) {
//            sb.append("TOTAL (+ zeros) outer-products = \n" + totalOp + "\n\n");
//            sb.append("TOTAL (- zeros) outer-products = \n" + totalNop + "\n\n");
//        }


        return sb.toString();
    }

    private static boolean allZero(double[] x) {
        boolean result = x[0] == 0.0;
        for (int i = 1; i < x.length && result; ++i) {
            result = x[i] == 0.0;
        }
        return result;
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

    public final ContinuousRateTransformation getRateTransformation() {
        return rateTransformation;
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

        final double[] tipPartial = forceCompletelyObserved ?
                dataModel.getTipPartial(tipIndex, true) :
                dataModel.getTipPartial(tipIndex);

//        if (precisionType == PrecisionType.SCALAR) {
//            System.err.println(new dr.math.matrixAlgebra.Vector(tipPartial));
//        }
//
//        final double[] tipPartial =
//                forceCompletelyObserved ?
//                dataModel.getTipPartial(tipIndex, true) :
//                        dataModel.getTipPartial(tipIndex);
//
//
//
//        // TODO Need specify the precision pattern for the returned partial
//
//        if (forceCompletelyObserved) {
//            tipPartial[dimTrait] = Double.POSITIVE_INFINITY;
////            System.err.println("FORCED");
////            System.exit(-1);
//        }
//
//        if (cdi instanceof ContinuousDiffusionIntegrator.Basic) {
//            System.err.println(tipPartial[dimTrait]);
//        }
//
////        System.err.println(cdi.getClass().getCanonicalName());

        setTipDataDirectly(tipIndex, tipPartial);
    }

    public void setTipDataDirectly(int tipIndex, double[] tipPartial) {
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

        int operationCount = nodeOperations.size();
        int k = 0;
        for (NodeOperation op : nodeOperations) {
            int nodeNum = op.getNodeNumber();

            operations[k + 0] = getActiveNodeIndex(op.getNodeNumber());
            operations[k + 1] = getActiveNodeIndex(op.getLeftChild());    // source node 1
            operations[k + 2] = getActiveMatrixIndex(op.getLeftChild());  // source matrix 1
            operations[k + 3] = getActiveNodeIndex(op.getRightChild());   // source node 2
            operations[k + 4] = getActiveMatrixIndex(op.getRightChild()); // source matrix 2

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

        cdi.updatePostOrderPartials(operations, operationCount, computeWishartStatistics);

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

    public final int getActiveNodeIndex(final int index) {
        return partialBufferHelper.getOffsetIndex(index);
    }

    public final int getActiveMatrixIndex(final int index) {
        return diffusionProcessDelegate.getMatrixIndex(index);
    }

    public void getPostOrderPartial(final int nodeNumber, double[] vector) {
        cdi.getPostOrderPartial(getActiveNodeIndex(nodeNumber), vector);
    }

    public void getPreOrderPartial(final int nodeNumber, double[] vector) {
        cdi.getPreOrderPartial(getActiveNodeIndex(nodeNumber), vector);
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
//        assert (callbackLikelihood != null);
//        callbackLikelihood.makeDirty();
//        computeWishartStatistics = true;
//        callbackLikelihood.getLogLikelihood();
//        computeWishartStatistics = false;
        return wishartStatistics;
    }

    public void setComputeWishartStatistics(boolean computeWishartStatistics) {
        this.computeWishartStatistics = computeWishartStatistics;
    }

    @Override
    public MatrixParameterInterface getPrecisionParamter() {
        return diffusionModel.getPrecisionParameter();
    }

    public ContinuousDataLikelihoodDelegate createObservedDataOnly(ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        return new ContinuousDataLikelihoodDelegate(likelihoodDelegate.tree,
                likelihoodDelegate.diffusionModel,
                likelihoodDelegate.dataModel,
                likelihoodDelegate.rootPrior,
                likelihoodDelegate.rateTransformation,
                likelihoodDelegate.rateModel,
                true);
    }

    private final static boolean DEBUG_OUTER_PRODUCTS = false;
    private final static boolean DEBUG_MISSING_DISTRIBUTION = true;

}
