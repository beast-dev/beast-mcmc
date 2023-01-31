/*
 * WishartStatisticsWrapper.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import dr.math.matrixAlgebra.Matrix;

import java.util.Arrays;
import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 */
public class WishartStatisticsWrapper extends AbstractModel implements ConjugateWishartStatisticsProvider, Loggable, Reportable {

    public static final String PARSER_NAME = "wishartStatistics";
    public static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    public WishartStatisticsWrapper(final String name,
                                    final String traitName,
                                    final TreeDataLikelihood dataLikelihood,
                                    final ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name);
        this.dataLikelihood = dataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;

        this.dimTrait = likelihoodDelegate.getTraitDim();
        this.numTrait = likelihoodDelegate.getTraitCount();
        this.tipCount = dataLikelihood.getTree().getExternalNodeCount();
        this.dimPartial = dimTrait + 1;

        dataLikelihood.addModelListener(this);
        dataLikelihood.addModel(this);
        String partialTraitName = getTipTraitName(traitName);
        tipSampleTrait = dataLikelihood.getTreeTrait(partialTraitName);

        treeTraversalDelegate = new LikelihoodTreeTraversal(
                dataLikelihood.getTree(),
                dataLikelihood.getBranchRateModel(),
                TreeTraversal.TraversalType.POST_ORDER);

        if (likelihoodDelegate.getIntegrator() instanceof MultivariateIntegrator) {

            ContinuousTraitPartialsProvider dataProvider = likelihoodDelegate.getDataModel();

            if (!dataProvider.suppliesWishartStatistics()) {
                dataProvider = new EmptyTraitDataModel(traitName, dataProvider.getParameter(),
                        dataProvider.getTraitDimension(), PrecisionType.SCALAR);
            }

            outerProductDelegate = ContinuousDataLikelihoodDelegate.createObservedDataOnly(
                    likelihoodDelegate, dataProvider);

        } else {
            outerProductDelegate = likelihoodDelegate;
        }

        traitDataKnown = false;
        outerProductsKnown = false;
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {

        if (!outerProductsKnown) {
            computeOuterProducts();
            outerProductsKnown = true;
        }

        return wishartStatistics;
    }

    public void simulateMissingTraits() {

        likelihoodDelegate.fireModelChanged(); // Force new sample!

        double[] sample = (double[]) tipSampleTrait.getTrait(dataLikelihood.getTree(), null);

        if (DEBUG) {
            System.err.println("Attempting to simulate missing traits");
            System.err.println(new Vector(sample));
        }

        final ContinuousDiffusionIntegrator cdi = outerProductDelegate.getIntegrator();
        assert (cdi instanceof ContinuousDiffusionIntegrator.Basic);

        double[] buffer = new double[dimPartial * numTrait];
        for (int trait = 0; trait < numTrait; ++trait) {
            buffer[trait * dimPartial + dimTrait] = Double.POSITIVE_INFINITY;
        }

        for (int tip = 0; tip < tipCount; ++tip) {
            int sampleOffset = tip * dimTrait * numTrait;
            int bufferOffset = 0;
            for (int trait = 0; trait < numTrait; ++trait) {
                System.arraycopy(sample, sampleOffset, buffer, bufferOffset, dimTrait);
                sampleOffset += dimTrait;
                bufferOffset += dimPartial;
            }
            outerProductDelegate.setTipDataDirectly(tip, buffer);
        }

        if (DEBUG) {
            System.err.println("Finished draw");
        }
    }

    private void computeOuterProducts() {

        // Make sure partials on tree are ready
        dataLikelihood.getLogLikelihood();

        if (likelihoodDelegate != outerProductDelegate) {
            simulateMissingTraits();
        }

        treeTraversalDelegate.updateAllNodes();
        treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

        final List<DataLikelihoodDelegate.BranchOperation> branchOperations = treeTraversalDelegate.getBranchOperations();
        final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = treeTraversalDelegate.getNodeOperations();

        final NodeRef root = dataLikelihood.getTree().getRoot();

        outerProductDelegate.setComputeWishartStatistics(true);
        outerProductDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());
        outerProductDelegate.setComputeWishartStatistics(false);

        wishartStatistics = outerProductDelegate.getWishartStatistics();

        if (DEBUG) {
            System.err.println("WS: " + wishartStatistics);
        }
    }

    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return likelihoodDelegate.getDiffusionModel().getPrecisionParameter();
    }

    @Override
    protected void storeState() {
        savedTraitDataKnown = traitDataKnown;
        savedOuterProductsKnown = outerProductsKnown;

        if (outerProductsKnown) {
            if (savedWishartStatistics == null) {
                savedWishartStatistics = wishartStatistics.clone();
            } else {
                wishartStatistics.copyTo(savedWishartStatistics);
            }
        }
    }

    @Override
    protected void restoreState() {
        traitDataKnown = savedTraitDataKnown;
        outerProductsKnown = savedOuterProductsKnown;

        if (outerProductsKnown) {
            WishartSufficientStatistics tmp = wishartStatistics;
            wishartStatistics = savedWishartStatistics;
            savedWishartStatistics = tmp;
        }
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        outerProductsKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        outerProductsKnown = false;
        // TODO If no partially missing traits and diffusion model hit, then no update necessary
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.hasId() ? xo.getId() : PARSER_NAME;
            String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

            final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

            if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("May not provide a sequence data likelihood in the precision Gibbs sampler");

            }

            final ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

            return new WishartStatisticsWrapper(name,
                    traitName,
                    treeDataLikelihood, continuousData);
        }

        /**
         * @return an array of syntax rules required by this element.
         * Order is not important.
         */
        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return syntax;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return WishartStatisticsWrapper.class;
        }

        /**
         * @return Parser name, which is identical to name of xml element parsed by it.
         */
        @Override
        public String getParserName() {
            return PARSER_NAME;
        }

        private final XMLSyntaxRule[] syntax = new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TRAIT_NAME, true),
        };
    };

    private final LikelihoodTreeTraversal treeTraversalDelegate;
    private final TreeTrait tipSampleTrait;

    private final int dimTrait;
    private final int numTrait;
    private final int tipCount;
    private final int dimPartial;

    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final ContinuousDataLikelihoodDelegate outerProductDelegate;
    private final TreeDataLikelihood dataLikelihood;

    private boolean traitDataKnown;
    private boolean outerProductsKnown;

    private boolean savedTraitDataKnown;
    private boolean savedOuterProductsKnown;

    private WishartSufficientStatistics wishartStatistics;
    private WishartSufficientStatistics savedWishartStatistics;

    private static final boolean DEBUG = false;

    @Override
    public LogColumn[] getColumns() {

        int sampleLength = 0;
        if (tipSampleTrait != null) {
            double[] sample = (double[]) tipSampleTrait.getTrait(dataLikelihood.getTree(), null);
            sampleLength = sample.length;
        }

        LogColumn[] columns = new LogColumn[dimTrait * dimTrait + sampleLength];

        int index = 0;
        for (int i = 0; i < dimTrait; ++i) {
            for (int j = 0; j < dimTrait; ++j) {
                columns[index] = new OuterProductColumn("OP" + (i + 1) + "" + (j + 1), index);
                ++index;
            }
        }

        for (int i = 0; i < sampleLength; ++i) {
            columns[index] = new TipSampleColumn("TIP" + (i + 1), i);
            ++index;
        }

        return columns;
    }

    @Override
    public String getReport() {

        WishartSufficientStatistics stats = getWishartStatistics();
        double[][] treePrec = likelihoodDelegate.getTreePrecision();


        int dimTrait = likelihoodDelegate.getTraitDim();
        int nTaxa = treePrec.length;

        double[] traits = (double[]) tipSampleTrait.getTrait(dataLikelihood.getTree(), null); // should just re-sample the existing traits

        DenseMatrix64F traitMat = DenseMatrix64F.wrap(nTaxa, dimTrait, traits);

        double[] rootMean = likelihoodDelegate.getRootPrior().getMean();
        DenseMatrix64F rootVec = DenseMatrix64F.wrap(dimTrait, 1, rootMean);
        double[] ones = new double[nTaxa];
        Arrays.fill(ones, 1.0);
        DenseMatrix64F oneVec = DenseMatrix64F.wrap(nTaxa, 1, ones);
        DenseMatrix64F rootMat = new DenseMatrix64F(nTaxa, dimTrait);

        CommonOps.multTransB(oneVec, rootVec, rootMat);
        CommonOps.add(traitMat, -1, rootMat);

        DenseMatrix64F treePrecMat = new DenseMatrix64F(nTaxa, nTaxa);
        for (int i = 0; i < nTaxa; i++) {
            for (int j = 0; j < nTaxa; j++) {
                treePrecMat.set(i, j, treePrec[i][j]);
            }
        }

        DenseMatrix64F bufferMat = new DenseMatrix64F(nTaxa, dimTrait);
        CommonOps.mult(treePrecMat, traitMat, bufferMat);

        DenseMatrix64F scaleMat = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransA(traitMat, bufferMat, scaleMat);

        StringBuilder sb = new StringBuilder("WishartStatisticsWrapper report:\n\n");

        sb.append("Scale matrix (naive):\n");
        sb.append(new Matrix(scaleMat.data, dimTrait, dimTrait));

        sb.append("\n");


        sb.append("Scale matrix (recursive):\n");
        sb.append(new Matrix(stats.getScaleMatrix(), dimTrait, dimTrait));

        sb.append("\n\n");

        return sb.toString();
    }

    private class OuterProductColumn extends NumberColumn {

        private int index;

        private OuterProductColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            WishartSufficientStatistics ws = getWishartStatistics();
            return ws.getScaleMatrix()[index];
        }
    }

    private class TipSampleColumn extends NumberColumn {

        private int index;

        private TipSampleColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            double[] sample = (double[]) tipSampleTrait.getTrait(dataLikelihood.getTree(), null);
            return sample[index];
        }
    }
}
