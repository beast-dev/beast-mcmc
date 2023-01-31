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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.ConditionalPrecisionAndTransform2;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.PermutationIndices;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */
public class BranchRateGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final int nTraits;
//    private final int dim;
    private final Parameter rateParameter;
    private final DifferentiableBranchRates branchRateModel;
    private final ContinuousTraitGradientForBranch branchProvider;

//    private final DenseMatrix64F matrix0;
//    private final DenseMatrix64F matrix1;
//    private final DenseMatrix64F vector0;

    public BranchRateGradient(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert (treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;

        BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        this.branchRateModel = (brm instanceof DifferentiableBranchRates) ? (DifferentiableBranchRates) brm : null;

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
        final int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
//
//        if (likelihoodDelegate.getIntegrator() instanceof SafeMultivariateWithDriftIntegrator) {
//            throw new RuntimeException("Not yet implemented with drift");
//        }

        branchProvider = new ContinuousTraitGradientForBranch.Default(dim);

//        matrix0 = new DenseMatrix64F(dim, dim);
//        matrix1 = new DenseMatrix64F(dim, dim);
//        vector0 = new DenseMatrix64F(dim, 1);
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

        treeDataLikelihood.makeDirty(); // TODO Remove after we figure out why this is necessary

        double[] result = new double[rateParameter.getDimension()];

        // TODO Do single call to traitProvider with node == null (get full tree)
//        List<BranchSufficientStatistics> statisticsForTree = (List<BranchSufficientStatistics>)
//                treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {

                List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);

                assert (statisticsForNode.size() == nTraits);

                final double rate = branchRateModel.getBranchRate(tree, node);
                final double differential = branchRateModel.getBranchRateDifferential(tree, node);
                final double scaling = differential / rate;

                double gradient = 0.0;
                for (int trait = 0; trait < nTraits; ++trait) {
                    gradient += branchProvider.getGradientForBranch(statisticsForNode.get(trait), scaling);
                }

                final int destinationIndex = getParameterIndexFromNode(node);
                assert (destinationIndex != -1);

                result[destinationIndex] = gradient;
            }
        }

        return result;
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

    //TODO: (Alex) Remove code duplication HERE:

    @Override
    public double[] getDiagonalHessianLogDensity() {

        double diagonalHessianResult[] = new double[rateParameter.getDimension()];

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
                assert (statisticsForNode.size() == nTraits);

                final double rate = branchRateModel.getBranchRate(tree, node);
                final double differential = branchRateModel.getBranchRateDifferential(tree, node);
                final double scaling = differential / rate;
                final double secondDifferential = branchRateModel.getBranchRateSecondDifferential(tree, node);
                final double secondScaling = secondDifferential / rate;

                double hessian = 0.0;

                for (int trait = 0; trait < nTraits; ++trait) {
                    hessian += getDiagonalHessianLogDensity(statisticsForNode.get(trait), scaling, secondScaling);
                }

                final int destinationIndex = getParameterIndexFromNode(node);
                assert (destinationIndex != -1);

                diagonalHessianResult[destinationIndex] = hessian;
            }
        }
        return diagonalHessianResult;
    }

    private double getDiagonalHessianLogDensity(BranchSufficientStatistics statistics,
                                                double differentialScaling,
                                                double secondDifferentialScaling) {
        final DenseMatrix64F matrix0;
        final DenseMatrix64F matrix1;
        final DenseMatrix64F matrix2;

        final DenseMatrix64F vector0;

        final int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();

        matrix0 = new DenseMatrix64F(dim, dim);
        matrix1 = new DenseMatrix64F(dim, dim);
        matrix2 = new DenseMatrix64F(dim, dim);

        vector0 = new DenseMatrix64F(dim, dim);

        final NormalSufficientStatistics below = statistics.getBelow();
        final NormalSufficientStatistics branch = statistics.getBranch();
        final NormalSufficientStatistics above = statistics.getAbove();

        NormalSufficientStatistics jointStatistics = ((ContinuousTraitGradientForBranch.Default) branchProvider).computeJointStatistics(below, above, dim);

        ((ContinuousTraitGradientForBranch.Default) branchProvider).makeDeltaVector(vector0, jointStatistics, above);
        DenseMatrix64F delta = vector0; // delta is expectation vector M_i - n_i

        DenseMatrix64F Qi = above.getRawPrecision();
        DenseMatrix64F Vi = jointStatistics.getRawVariance();

        ((ContinuousTraitGradientForBranch.Default) branchProvider).makeGradientMatrices0(matrix1, matrix0, statistics, differentialScaling);
        DenseMatrix64F logDetComponent = matrix0;
        DenseMatrix64F quadraticComponent = matrix1;

        DenseMatrix64F squareLogDetComponent = matrix2;
        CommonOps.mult(logDetComponent, logDetComponent, squareLogDetComponent);

        //The diagonal hessian for observed data only requires hess0 and
        double hess0 = 0.0;
        for (int row = 0; row < dim; ++row) {
            hess0 += 0.5 * squareLogDetComponent.unsafe_get(row, row);
        }

        DenseMatrix64F dQuadraticComponent1 = matrix0;
        CommonOps.mult(squareLogDetComponent, Qi, dQuadraticComponent1);

        DenseMatrix64F VdQuadraticComponent1 = matrix2;
        CommonOps.mult(Vi, dQuadraticComponent1, VdQuadraticComponent1);

        double hess1 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                hess1 -= delta.unsafe_get(row, 0)
                        * dQuadraticComponent1.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);

            }
            hess1 -= VdQuadraticComponent1.unsafe_get(row, row);
        }

        ((ContinuousTraitGradientForBranch.Default) branchProvider).makeGradientMatrices1(matrix2, quadraticComponent, jointStatistics);
        DenseMatrix64F additionalVariance = matrix2;

        DenseMatrix64F quadrupleComponent = matrix0;
        CommonOps.mult(additionalVariance, quadraticComponent, quadrupleComponent); //ups V ups

        DenseMatrix64F quadraticVariance = matrix2;
        CommonOps.mult(quadrupleComponent, Vi, quadraticVariance); //ups V ups V

        double hess2 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                hess2 += delta.unsafe_get(row, 0)
                        * quadrupleComponent.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);

            }
            hess2 += 0.5 * quadraticVariance.unsafe_get(row, row);
        }

        DenseMatrix64F Si2 = matrix0;
        CommonOps.scale(secondDifferentialScaling, branch.getRawVariance(), Si2);

        DenseMatrix64F secondLogDetComponent = matrix2;
        CommonOps.mult(Qi, Si2, secondLogDetComponent);

        DenseMatrix64F secondQuadraticComponent = matrix0;
        CommonOps.mult(secondLogDetComponent, Qi, secondQuadraticComponent); // Q Si2 Q

        DenseMatrix64F VsecondQuadraticComponent = matrix1;
        CommonOps.mult(Vi, secondQuadraticComponent, VsecondQuadraticComponent);

        double hess3 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                hess3 += delta.unsafe_get(row, 0)
                        * secondQuadraticComponent.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);
            }
            hess3 -= 0.5 * secondLogDetComponent.unsafe_get(row, row);
            hess3 += 0.5 * VsecondQuadraticComponent.unsafe_get(row, row);
        }

        return hess0 + hess1 + hess2 + hess3;
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }

    interface ContinuousTraitGradientForBranch {

        double getGradientForBranch(BranchSufficientStatistics statistics, double differentialScaling);

        class Default implements ContinuousTraitGradientForBranch {

            private final DenseMatrix64F matrix0;
            private final DenseMatrix64F matrix1;
            private final DenseMatrix64F vector0;

            private final int dim;

            public Default(int dim) {
                this.dim = dim;

                matrix0 = new DenseMatrix64F(dim, dim);
                matrix1 = new DenseMatrix64F(dim, dim);
                vector0 = new DenseMatrix64F(dim, 1);
            }

            @Override
            public double getGradientForBranch(BranchSufficientStatistics statistics, double differentialScaling) {

                final NormalSufficientStatistics below = statistics.getBelow();
                final NormalSufficientStatistics branch = statistics.getBranch();
                final NormalSufficientStatistics above = statistics.getAbove();

                if (DEBUG) {
                    System.err.println("B = " + statistics.toVectorizedString());
                }

//                if (DEBUG) {
//                    System.err.println("\tQi = " + NormalSufficientStatistics.toVectorizedString(Qi));
//                    System.err.println("\tSi = " + NormalSufficientStatistics.toVectorizedString(Si));
//                }

                DenseMatrix64F Qi = above.getRawPrecision();

                DenseMatrix64F logDetComponent = matrix0;
                DenseMatrix64F quadraticComponent = matrix1;
                makeGradientMatrices0(quadraticComponent, logDetComponent, statistics, differentialScaling);

                double grad1 = 0.0;
                for (int row = 0; row < dim; ++row) {
                    grad1 -= 0.5 * logDetComponent.unsafe_get(row, row);
                }

                if (DEBUG) {
                    System.err.println("grad1 = " + grad1);
                }


//                if (DEBUG) {
//                    System.err.println("\tQi = " + NormalSufficientStatistics.toVectorizedString(Qi));
//                    System.err.println("\tQQ = " + NormalSufficientStatistics.toVectorizedString(quadraticComponent));
//                }

                NormalSufficientStatistics jointStatistics = computeJointStatistics(below, above, dim);

                DenseMatrix64F additionalVariance = matrix0; //new DenseMatrix64F(dim, dim);
                makeGradientMatrices1(additionalVariance, quadraticComponent,
                        jointStatistics);


                DenseMatrix64F delta = vector0;
                makeDeltaVector(delta, jointStatistics, above);

                double grad2 = 0.0;
                for (int row = 0; row < dim; ++row) {
                    for (int col = 0; col < dim; ++col) {

                        grad2 += 0.5 * delta.unsafe_get(row, 0)
                                * quadraticComponent.unsafe_get(row, col)
                                * delta.unsafe_get(col, 0);

                    }

                    grad2 += 0.5 * additionalVariance.unsafe_get(row, row);
                }

                if (DEBUG) {
                    System.err.println("\tjoint mean  = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawMean()));
                    System.err.println("\tabove mean = " + NormalSufficientStatistics.toVectorizedString(above.getRawMean()));
                    System.err.println("\tbelow mean  = " + NormalSufficientStatistics.toVectorizedString(below.getRawMean()));
                    System.err.println("\tjoint precision  = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawPrecision()));
                    System.err.println("\tabove precision = " + NormalSufficientStatistics.toVectorizedString(above.getRawPrecision()));
                    System.err.println("\tbelow precision  = " + NormalSufficientStatistics.toVectorizedString(below.getRawPrecision()));
                    System.err.println("\tbelow variance   = " + NormalSufficientStatistics.toVectorizedString(below.getRawVariance()));
                    System.err.println("\tquadratic      = " + NormalSufficientStatistics.toVectorizedString(quadraticComponent));
                    System.err.println("\tadditional     = " + NormalSufficientStatistics.toVectorizedString(additionalVariance));
                    System.err.println("delta: " + NormalSufficientStatistics.toVectorizedString(delta));
                    System.err.println("grad2 = " + grad2);
                }

                // W.r.t. drift
                // TODO: Fix delegate to (possibly) un-link drift from arbitrary rate
                DenseMatrix64F Di = new DenseMatrix64F(dim, 1);
                CommonOps.scale(differentialScaling, branch.getRawMean(), Di);

                double grad3 = 0.0;
                for (int row = 0; row < dim; ++row) {
                    for (int col = 0; col < dim; ++col) {

                        grad3 += delta.unsafe_get(row, 0)
                                * Qi.unsafe_get(row, col)
                                * Di.unsafe_get(col, 0);

                    }
                }

                if (DEBUG) {
                    System.err.println("\tDi     = " + NormalSufficientStatistics.toVectorizedString(branch.getRawMean()));
                    System.err.println("grad3 = " + grad3);
                }

                return grad1 + grad2 + grad3;

            }

            public static NormalSufficientStatistics computeJointStatistics(NormalSufficientStatistics below,
                                                                            NormalSufficientStatistics above,
                                                                            int dim) {

                PermutationIndices indices = new PermutationIndices(below.getRawPrecision());

                if (indices.getNumberOfInfiniteDiagonals() == dim) {

                    return computeJointFullyObserved(below, dim);

                } else if (indices.getNumberOfZeroDiagonals() == dim) {

                    return computeJointFullyMissing(above, dim);

                } else if (indices.getNumberOfZeroDiagonals() == 0 || indices.getNumberOfInfiniteDiagonals() == 0) {

                    return computeJointLatent(below, above, dim);

                } else {

                    return computeJointPartiallyMissing(below, above, indices, dim);
                }
            }

            private static NormalSufficientStatistics computeJointFullyObserved(NormalSufficientStatistics below,
                                                                                int dim) {

                return new NormalSufficientStatistics(
                        below.getRawMean(), below.getRawPrecision(), new DenseMatrix64F(dim, dim));
            }

            private static NormalSufficientStatistics computeJointFullyMissing(NormalSufficientStatistics above,
                                                                               int dim) {

                return new NormalSufficientStatistics(
                        above.getRawMean(), above.getRawPrecision(), above.getRawVariance());
            }

            private static NormalSufficientStatistics computeJointLatent(NormalSufficientStatistics below,
                                                                         NormalSufficientStatistics above,
                                                                         int dim) {

                DenseMatrix64F mean = new DenseMatrix64F(dim, 1);
                DenseMatrix64F precision = new DenseMatrix64F(dim, dim);
                DenseMatrix64F variance = new DenseMatrix64F(dim, dim);

                CommonOps.add(below.getRawPrecision(), above.getRawPrecision(), precision);
                safeInvert2(precision, variance, false);

                safeWeightedAverage(
                        new WrappedVector.Raw(below.getRawMean().getData(), 0, dim),
                        below.getRawPrecision(),
                        new WrappedVector.Raw(above.getRawMean().getData(), 0, dim),
                        above.getRawPrecision(),
                        new WrappedVector.Raw(mean.getData(), 0, dim),
                        variance,
                        dim);

                return new NormalSufficientStatistics(mean, precision, variance);
            }

            private static NormalSufficientStatistics computeJointPartiallyMissing(NormalSufficientStatistics below,
                                                                                   NormalSufficientStatistics above,
                                                                                   PermutationIndices indices,
                                                                                   int dim) {

                DenseMatrix64F mean = new DenseMatrix64F(dim, 1);
                DenseMatrix64F precision = new DenseMatrix64F(dim, dim);
                DenseMatrix64F variance = new DenseMatrix64F(dim, dim);

                if (indices.getNumberOfNonZeroFiniteDiagonals() != 0) {
                    throw new RuntimeException("Unsure if this works for latent trait below");
                    // TODO Probably need to add in below.precision somewhere below
                }

                ConditionalPrecisionAndTransform2 transform = new ConditionalPrecisionAndTransform2(
                        above.getRawPrecision(),
                        indices.getZeroIndices(),
                        indices.getInfiniteIndices()
                );

                double[] result = transform.getConditionalMean(below.getRawMean().getData(), 0,
                        above.getRawMean().getData(), 0);

                copyRowsAndColumns(above.getRawPrecision(), precision,
                        indices.getZeroIndices(), indices.getZeroIndices(), false);

                scatterRowsAndColumns(transform.getConditionalVariance(), variance,
                        indices.getZeroIndices(), indices.getZeroIndices(), false);

                int index = 0;
                for (int zero : indices.getZeroIndices()) {
                    mean.unsafe_set(zero, 0, result[index++]);
                }

                for (int infinite : indices.getInfiniteIndices()) {
                    mean.unsafe_set(infinite, 0, below.getMean(infinite));
                    precision.unsafe_set(infinite, infinite, Double.POSITIVE_INFINITY);
                }

                return new NormalSufficientStatistics(mean, precision, variance);
            }

            public void makeGradientMatrices0(DenseMatrix64F matrix1, DenseMatrix64F logDetComponent,
                                              BranchSufficientStatistics statistics, double differentialScaling) {

                final NormalSufficientStatistics above = statistics.getAbove();
                final NormalSufficientStatistics branch = statistics.getBranch();

                DenseMatrix64F Qi = above.getRawPrecision();
                CommonOps.scale(differentialScaling, branch.getRawVariance(), matrix1); //matrix1 = Si
                CommonOps.mult(Qi, matrix1, logDetComponent); //matrix0 = logDetComponent
                CommonOps.mult(logDetComponent, Qi, matrix1); //matrix1 = QuadraticComponent

            }

            public void makeGradientMatrices1(DenseMatrix64F additionalVariance, DenseMatrix64F quadraticComponent,
                                              NormalSufficientStatistics jointStatistics) {

                CommonOps.mult(quadraticComponent, jointStatistics.getRawVariance(), additionalVariance);
            }

            public void makeDeltaVector(DenseMatrix64F delta, NormalSufficientStatistics jointStatistics,
                                        NormalSufficientStatistics above) {
                for (int row = 0; row < dim; ++row) {
                    delta.unsafe_set(row, 0,
                            jointStatistics.getRawMean().unsafe_get(row, 0) - above.getMean(row)
                    );
                }
            }

        }
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

    public double[] getNumericalGradient() {
        double[] savedValues = rateParameter.getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, rateParameter.getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            rateParameter.setParameterValue(i, savedValues[i]);
        }

        return testGradient;
    }

    @Override
    public String getReport() {

        double[] testGradient = getNumericalGradient();

        NumericalHessianFromGradient numericalHessianFromGradient = new NumericalHessianFromGradient(this);

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

        sb.append("Peeling diagonal hessian: ").append(new dr.math.matrixAlgebra.Vector(getDiagonalHessianLogDensity()));
        sb.append("\n");
        sb.append("numeric diagonal hessian: ").append(new dr.math.matrixAlgebra.Vector(NumericalDerivative.diagonalHessian(numeric1, getParameter().getParameterValues())));
        sb.append("\n");
        sb.append("Another numeric diagonal hessian: ").append(new dr.math.matrixAlgebra.Vector(numericalHessianFromGradient.getDiagonalHessianLogDensity()));
        sb.append("\n");

        return sb.toString();
    }

    private static final boolean DEBUG = false;

    @Override
    public LogColumn[] getColumns() {

        return Loggable.getColumnsFromReport(this, "BranchRateGradient report");
    }
}
