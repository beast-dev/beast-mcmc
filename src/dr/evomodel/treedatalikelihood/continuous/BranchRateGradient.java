/*
 * TreeTipGradient.java
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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.UnivariateFunction;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BranchRateGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final int nTraits;
    private final int dim;
    private final Parameter rateParameter;
    private final ArbitraryBranchRates branchRateModel;
    private final MatrixParameterInterface diffusionParameter;

    private final DenseMatrix64F L;
    private final DenseMatrix64F Lt;
    private final DenseMatrix64F matrix0;
    private final DenseMatrix64F matrix1;
    private final DenseMatrix64F vector0;

    public BranchRateGradient(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert(treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;
        this.diffusionParameter = likelihoodDelegate.getDiffusionModel().getPrecisionParameter();

        BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        this.branchRateModel = (brm instanceof ArbitraryBranchRates) ? (ArbitraryBranchRates) brm : null;

        // TODO Move into different constructor / parser
        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        treeTraitProvider = unchecked;

        assert (treeTraitProvider != null);

        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }
        dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        assert (dim == diffusionParameter.getColumnDimension());

        if (likelihoodDelegate.getIntegrator() instanceof SafeMultivariateWithDriftIntegrator) {
            throw new RuntimeException("Not yet implemented with drift");
        }

        L = new DenseMatrix64F(dim, dim);
        Lt = new DenseMatrix64F(dim, dim);
        matrix0 = new DenseMatrix64F(dim, dim);
        matrix1 = new DenseMatrix64F(dim, dim);
        vector0 = new DenseMatrix64F(dim, 1);
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        getCholeskyDecomposition(diffusionParameter, L, Lt);

        // TODO Do single call to traitProvider with node == null (get full tree)
//        List<BranchSufficientStatistics> statisticsForTree = (List<BranchSufficientStatistics>)
//                treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {

                List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);

                assert (statisticsForNode.size() == nTraits);

                double gradient = 0.0;

//                double rate = branchRateModel.getBranchRate(tree, node);
                double differential = branchRateModel.getBranchRateDifferential(tree, node);

                for (int trait = 0; trait < nTraits; ++trait) {
                    gradient += //getGradientForBranchPrecision(rate)
//                            +
                                    getGradientForBranchExponential(statisticsForNode.get(trait), L, Lt, differential);
                    // TODO Fix both above for differ rate parametrization
                }

                final int destinationIndex = getParameterIndexFromNode(node);
                assert (destinationIndex != -1);

                result[destinationIndex] = gradient;
            }
        }

        return result;
    }

    private void getCholeskyDecomposition(MatrixParameterInterface input, DenseMatrix64F L, DenseMatrix64F Lt) {

        CholeskyDecomposition<DenseMatrix64F> cholesky = DecompositionFactory.chol(dim, true);
        DenseMatrix64F matrix = DenseMatrix64F.wrap(dim, dim, input.getParameterValues());

        if( !cholesky.decompose(matrix)) {
            throw new RuntimeException("Cholesky decomposition failed!");
        }

//        return cholesky.getT(null);
        cholesky.getT(L);
        CommonOps.transpose(L, Lt);
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

    private double getGradientForBranchPrecision(double rate) {
        return  0.0; // 0.5 / rate;
    }

    private double getGradientForBranchExponential(BranchSufficientStatistics statistics,
                                                   DenseMatrix64F L, DenseMatrix64F Lt,
                                                   double differentialScaling) {

        final NormalSufficientStatistics child = statistics.getChild();
        final NormalSufficientStatistics branch = statistics.getBranch();
        final NormalSufficientStatistics parent = statistics.getParent();

        if (DEBUG) {
//            System.err.println("L = " + NormalSufficientStatistics.toVectorizedString(L));
            System.err.println("B = " + statistics.toVectorizedString());
        }

        DenseMatrix64F Qi = parent.getRawPrecision();
        DenseMatrix64F Si = matrix1;
        CommonOps.scale(differentialScaling, branch.getRawVariance(), Si);

        DenseMatrix64F logDetComponent = matrix0;
        CommonOps.mult(Qi, Si, logDetComponent);

        double grad1 = 0.0;
        for (int row = 0; row < dim; ++row) {
            grad1 -= 0.5 * logDetComponent.unsafe_get(row, row);
        }

        if (DEBUG) {
            System.err.println("grad1 = " + grad1);
        }

        DenseMatrix64F quadraticComponent = matrix1;
        CommonOps.mult(logDetComponent, Qi, quadraticComponent);

        DenseMatrix64F delta = vector0;

        for (int row = 0; row < dim; ++row) {
            delta.unsafe_set(row, 0,
                    child.getMean(row) - parent.getMean(row) // TODO Correct for drift
            );
        }

        double grad2 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                grad2 += 0.5 * delta.unsafe_get(row, 0)
                        * quadraticComponent.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);
//
//                V.unsafe_set(row, col,
//                        child.getVariance(row, col)
//                                - 2 * branch.getVariance(row, col)
//                                + parent.getVariance(row, col)
//                                + delta.unsafe_get(row, 0) * delta.unsafe_get(col, 0)
//                );
            }
        }

        if (DEBUG) {
            System.err.println("delta: " + NormalSufficientStatistics.toVectorizedString(delta));
            System.err.println("grad2 = " + grad2);
        }
//
//        if (DEBUG) {
//            System.err.println("V = " + NormalSufficientStatistics.toVectorizedString(V));
//        }
//
//        DenseMatrix64F tmp = matrix1;
//        CommonOps.mult(Lt, V, tmp);
//        DenseMatrix64F LtVL = matrix0;
//        CommonOps.mult(tmp, L, LtVL);
//
//        double gradient = 0.0;
//        for (int row = 0; row < dim; ++row) {
//            gradient += LtVL.unsafe_get(row, row);
//        }
//
//        if (DEBUG) {
//            System.err.println("LtVL = " + NormalSufficientStatistics.toVectorizedString(LtVL));
//        }

        if (DEBUG) {
            System.err.println("return: " + (grad1 + grad2));
        }

        return grad1 + grad2;
    }



    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                rateParameter.setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return rateParameter.getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return 0;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
        }
    };

//    private UnivariateFunction func = new UnivariateFunction(1) {
//
//
//
//        private int index;
//
//        void setIndex(int index) {
//            this.index = index;
//        }
//
//        @Override
//        public double evaluate(double argument) {
//
//            rateParameter.setParameterValue(index, argument);
//
////            NodeRef node = tree.getNode(branchRateModel.getNodeNumberFromParameterIndex(index));
////            List<BranchSufficientStatistics> statistics = treeTraitProvider.getTrait(tree, node);
////
////            assert (statistics.size() == 1);
//
////            MultivariateNormalDistribution child = new MultivariateNormalDistribution()
//
//
//            return treeDataLikelihood.getLogLikelihood();
//        }
//
//        @Override
//        public double getLowerBound() {
//            return 0;
//        }
//
//        @Override
//        public double getUpperBound() {
//            return Double.POSITIVE_INFINITY;
//        }
//    };


    @Override
    public String getReport() {



        double[] testGradient = NumericalDerivative.gradient(numeric1, rateParameter.getParameterValues());

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: " + (new dr.math.matrixAlgebra.Vector(getGradientLogDensity())));
        sb.append("\n");
        sb.append("numeric: " + (new dr.math.matrixAlgebra.Vector(testGradient)));
        sb.append("\n");

        final int index = 0;

        UnivariateFunction func = new UnivariateFunction() {


                @Override
                public double evaluate(double argument) {

                    rateParameter.setParameterValue(index, argument);

//                    NodeRef node = tree.getNode(branchRateModel.getNodeNumberFromParameterIndex(index));
//                    List<BranchSufficientStatistics> statistics = treeTraitProvider.getTrait(tree, node);
//
//                    NormalSufficientStatistics child = statistics.get(0).getChild();
//                    NormalSufficientStatistics branch = statistics.get(0).getBranch();
//                    NormalSufficientStatistics parent = statistics.get(0).getParent();
//
//                    DenseMatrix64F variance = new DenseMatrix64F(2 * dim, 2 * dim);
//
//                    for (int row = 0; row < dim; ++row) {
//                        for (int col = 0; col < dim; ++col) {
//                            variance.unsafe_set(row, col, child.getVariance(row, col));
//                            variance.unsafe_set(dim + row, dim + col, parent.getVariance(row, col));
//                            variance.unsafe_set(dim + row, col, branch.getVariance(row, col));
//                            variance.unsafe_set(row, dim + col, branch.getVariance(row, col));
//                        }
//                    }
//
//                    System.err.println(statistics.get(0).toVectorizedString());
//                    System.err.println(variance);
//                    System.exit(-1);
        //
        //            assert (statistics.size() == 1);

        //            MultivariateNormalDistribution child = new MultivariateNormalDistribution()

                    double precision = 2;
                    double sd = Math.sqrt((1 + 1.0 / argument) / precision);

                    NormalDistribution normal = new NormalDistribution(2, sd);

//                    treeDataLikelihood.makeDirty();
//                    return treeDataLikelihood.getLogLikelihood();
                    return normal.logPdf(10);
                }

                @Override
                public double getLowerBound() {
                    return 0;
                }

                @Override
                public double getUpperBound() {
                    return Double.POSITIVE_INFINITY;
                }
            };


        System.err.println(NumericalDerivative.firstDerivative(func, 3.0));

        return sb.toString();
    }

    private static final boolean DEBUG = true;
}
