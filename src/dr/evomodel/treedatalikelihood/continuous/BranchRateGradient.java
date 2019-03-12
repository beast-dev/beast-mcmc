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
import dr.evomodel.substmodel.DifferentiableSubstitutionModelUtil;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;
import no.uib.cipr.matrix.DenseMatrix;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert;
import static dr.math.matrixAlgebra.missingData.MissingOps.safeWeightedAverage;

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
    private final ArbitraryBranchRates branchRateModel;
    private final ContinuousTraitGradientForBranch branchProvider;

//    private final DenseMatrix64F matrix0;
//    private final DenseMatrix64F matrix1;
//    private final DenseMatrix64F vector0;

    public BranchRateGradient(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert(treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;

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
// TODO: (Alex) remove major code duplication of diag hessian
    private double getDiagonalHessianLogDensity(BranchSufficientStatistics statistics,
                                                double differentialScaling,
                                                double secondDifferentialScaling){
        final DenseMatrix64F matrix0;
        final DenseMatrix64F matrix1;
        final DenseMatrix64F matrix2;
        final DenseMatrix64F squareLogDetComponent;
        final DenseMatrix64F matrix3;

        final DenseMatrix64F vector0;


        final int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();

        matrix0 = new DenseMatrix64F(dim, dim);
        matrix1 = new DenseMatrix64F(dim, dim);
        matrix2 = new DenseMatrix64F(dim, dim);
        squareLogDetComponent = new DenseMatrix64F(dim, dim);
        matrix3 = new DenseMatrix64F(dim, dim);
        vector0 = new DenseMatrix64F(dim, dim);

        final NormalSufficientStatistics child = statistics.getChild();
        final NormalSufficientStatistics branch = statistics.getBranch();
        final NormalSufficientStatistics parent = statistics.getParent();

        DenseMatrix64F Qi = parent.getRawPrecision();
        DenseMatrix64F Si = matrix1; //new DenseMatrix64F(dim, dim); // matrix1;
        CommonOps.scale(differentialScaling, branch.getRawVariance(), Si);

        DenseMatrix64F logDetComponent = matrix0; //new DenseMatrix64F(dim, dim); //matrix0;
        CommonOps.mult(Qi, Si, logDetComponent);

        DenseMatrix64F quadraticComponent = matrix2; //new DenseMatrix64F(dim, dim); //matrix1;
        CommonOps.mult(logDetComponent, Qi, quadraticComponent);

        NormalSufficientStatistics jointStatistics = ((ContinuousTraitGradientForBranch.Default) branchProvider).computeJointStatistics(child, parent);

        // Mi - ni vector = delta
        DenseMatrix64F delta = vector0;
        for (int row = 0; row < dim; ++row) {
            delta.unsafe_set(row, 0,
                    jointStatistics.getRawMean().unsafe_get(row, 0) - parent.getMean(row)
            );
        }

        //New to Hessian Section
        //dUps
//        CommonOps.mult(quadraticComponent, Si, dUps);

         // TODO: dUps should also have a Qd^2 Sigma Q term (0 in this model.)
        DenseMatrix64F dUps = matrix3;
        CommonOps.scale(-2.0 * differentialScaling, branch.getRawVariance(), dUps);
        DifferentiableSubstitutionModelUtil.getTripleMatrixMultiplication(dim,
                new WrappedMatrix.WrappedDenseMatrix(quadraticComponent),
                new WrappedMatrix.WrappedDenseMatrix(dUps),
                new WrappedMatrix.WrappedDenseMatrix(Qi));

        CommonOps.mult(quadraticComponent, Si, squareLogDetComponent);


        //tip hessian
        double hess1 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                hess1 += 0.5 * delta.unsafe_get(row, 0)
                        * dUps.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);

            }
            hess1 += 0.5 * squareLogDetComponent.unsafe_get(row, row);
        }

        // don't need Si anymore since it only comes into play in both squareLogDetComonent and quadraticComponent

        DenseMatrix64F VdUps = matrix1;
        CommonOps.mult(jointStatistics.getRawVariance(), dUps, VdUps);


        double hess2 = 0.0;
        for (int row = 0; row < dim; ++row) {
            hess2 += 0.5 * VdUps.unsafe_get(row, row);
        }
        // now don't need dUps or VdUps

        DenseMatrix64F Si2 = matrix1;
        CommonOps.scale(secondDifferentialScaling ,branch.getRawVariance(), Si2);

        DenseMatrix64F secondDerivativeComponent = matrix3;
        CommonOps.mult(Qi, Si2, secondDerivativeComponent);

        double hess3 = 0.0;
        for (int row = 0; row < dim; ++row) {
            hess3 -= 0.5 * secondDerivativeComponent.unsafe_get(row, row);
        }

        // last place Qi will be used:
        DenseMatrix64F doublyQuadratic = matrix1;

        CommonOps.scale(1, jointStatistics.getRawVariance(), doublyQuadratic);
        DifferentiableSubstitutionModelUtil.getTripleMatrixMultiplication(dim,
                new WrappedMatrix.WrappedDenseMatrix(quadraticComponent),
                new WrappedMatrix.WrappedDenseMatrix(doublyQuadratic),
                new WrappedMatrix.WrappedDenseMatrix(quadraticComponent));

        DenseMatrix64F quadraticVariance = matrix3;
        CommonOps.mult(doublyQuadratic, jointStatistics.getRawVariance(), quadraticVariance);

        double hess4 = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                hess4 +=  delta.unsafe_get(row, 0)
                        * doublyQuadratic.unsafe_get(row, col)
                        * delta.unsafe_get(col, 0);

            }
            hess4 += 0.5 * quadraticVariance.unsafe_get(row, row);
        }




        return hess1 + hess2 + hess3 + hess4;
    }



//    private double getDiagonalHessianLogDensity(BranchSufficientStatistics statistics,
//                                                  double differentialScaling,
//                                                  double secondDifferentialScaling){
////                                                ContinuousTraitGradientForBranch.Default continuousTraitGradientForBranch) {
////                                                ContinuousTraitGradientForBranch continuousTraitGradientForBranch) {
//
//        final DenseMatrix64F matrix0;
//        final DenseMatrix64F matrix1;
//        final DenseMatrix64F matrix2;
//        final DenseMatrix64F matrix3;
//        final DenseMatrix64F vector0;
//        final DenseMatrix64F matrix99;
//        final DenseMatrix64F matrix4;
//        final DenseMatrix64F matrix5;
//        final DenseMatrix64F matrix6;
//        final DenseMatrix64F matrix7;
//        final DenseMatrix64F matrix8;
//        final DenseMatrix64F matrix9;
//        final DenseMatrix64F matrix10;
//        final DenseMatrix64F matrix20;
//        final DenseMatrix64F matrix21;
//        final int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
//
//        matrix0 = new DenseMatrix64F(dim, dim);
//        matrix1 = new DenseMatrix64F(dim, dim);
//        matrix2 = new DenseMatrix64F(dim, dim);
//        matrix3 = new DenseMatrix64F(dim, dim);
//        matrix99 = new DenseMatrix64F(dim, dim);
//        matrix4 = new DenseMatrix64F(dim, dim);
//        matrix5 = new DenseMatrix64F(dim, dim);
//        matrix6 = new DenseMatrix64F(dim, dim);
//        matrix7 = new DenseMatrix64F(dim, dim);
//        matrix8 = new DenseMatrix64F(dim, dim);
//        matrix9 = new DenseMatrix64F(dim, dim);
//        matrix10 = new DenseMatrix64F(dim, dim);
//        matrix20 = new DenseMatrix64F(dim, dim);
//        matrix21 = new DenseMatrix64F(dim, dim);
//        vector0 = new DenseMatrix64F(dim, 1);
//
//        final NormalSufficientStatistics child = statistics.getChild();
//        final NormalSufficientStatistics branch = statistics.getBranch();
//        final NormalSufficientStatistics parent = statistics.getParent();
//
//        DenseMatrix64F Qi = parent.getRawPrecision();
//        DenseMatrix64F Si = matrix1; //new DenseMatrix64F(dim, dim); // matrix1;
//        CommonOps.scale(differentialScaling, branch.getRawVariance(), Si);
//
//        DenseMatrix64F logDetComponent = matrix0; //new DenseMatrix64F(dim, dim); //matrix0;
//        CommonOps.mult(Qi, Si, logDetComponent);
//
//        DenseMatrix64F quadraticComponent = matrix1; //new DenseMatrix64F(dim, dim); //matrix1;
//        CommonOps.mult(logDetComponent, Qi, quadraticComponent);
//
//        NormalSufficientStatistics jointStatistics = ((ContinuousTraitGradientForBranch.Default) branchProvider).computeJointStatistics(child, parent);
//
//        DenseMatrix64F additionalVariance = matrix99; //new DenseMatrix64F(dim, dim);
//        CommonOps.mult(quadraticComponent, jointStatistics.getRawVariance(), additionalVariance);
//
//        DenseMatrix64F delta = vector0;
//        for (int row = 0; row < dim; ++row) {
//            delta.unsafe_set(row, 0,
//                    jointStatistics.getRawMean().unsafe_get(row, 0) - parent.getMean(row)
//            );
//        }
//
//        // new to the hessian section
//
//        DenseMatrix64F additionalQuadratic = matrix2;
//        CommonOps.mult(quadraticComponent, additionalVariance, additionalQuadratic);
//
//        DenseMatrix64F quadraticVariance = matrix3;
//        CommonOps.mult(additionalQuadratic, jointStatistics.getRawVariance(), quadraticVariance);
//
//        DenseMatrix64F R1 = matrix4;
//        CommonOps.mult(quadraticComponent, Si, R1);
//
//
//        DenseMatrix64F Rtemp = matrix5;
//        CommonOps.scale(secondDifferentialScaling ,branch.getRawVariance(), Rtemp);
//
//        DenseMatrix64F R2 = matrix6; // Qi * d^2 Sigma
//        CommonOps.mult(Qi, Rtemp, R2);
//
//        DenseMatrix64F BigA1 = matrix7;
//        CommonOps.mult(R1 ,Qi, BigA1);
//        CommonOps.scale(-2.0, BigA1, BigA1);
//
//        DenseMatrix64F A2 = matrix8;
//        CommonOps.mult(R2, Qi, A2);
//
//        DenseMatrix64F BigA = matrix10;
//        CommonOps.add(BigA1, A2, BigA);
//
//        DenseMatrix64F VdY = matrix9;
//        CommonOps.mult(jointStatistics.getRawVariance(), BigA, VdY);
//
//        DenseMatrix64F preAlternative = matrix20;
//        CommonOps.mult(logDetComponent, logDetComponent, preAlternative);
//
//        DenseMatrix64F Alternative = matrix21;
//        CommonOps.mult(preAlternative, Qi, Alternative);
//
//        double hess1 = 0.0;
//        for (int row = 0; row < dim; ++row) {
//            hess1 += 0.5 * quadraticVariance.unsafe_get(row, row);
//        }
//
//
//        double hess2 = 0.0;
//        for (int row = 0; row < dim; ++row) {
//            for (int col = 0; col < dim; ++col) {
//
//                hess2 += delta.unsafe_get(row, 0)
//                        * additionalQuadratic.unsafe_get(row, col)
//                        * delta.unsafe_get(col, 0);
//            }
//        }
//
//
//        double hess3 = 0.0; // (2)
//                for (int row = 0; row < dim; ++row) {
//                    hess3 += 0.5 * R1.unsafe_get(row, row);
//                }
//        double hess4 = 0.0; // (3)
//                for (int row = 0; row < dim; ++row) {
//                    hess4 -= 0.5 * R2.unsafe_get(row, row);
//                }
//
//        double hess5 = 0.0; // ( 1 )
//                for (int row = 0; row < dim; ++row) {
//                    for (int col = 0; col < dim; ++col) {
//                        hess5 += 0.5 * delta.unsafe_get(row, 0)
//                                * BigA.unsafe_get(row, col)
//                                * delta.unsafe_get(col, 0);
//                    }
//                    hess5 += 0.5 * VdY.unsafe_get(row, row);
//                }
//
//            double hess = hess1 + hess2 + hess3 + hess4 + hess5;
//            return hess;
//    }

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

                final NormalSufficientStatistics child = statistics.getChild();
                final NormalSufficientStatistics branch = statistics.getBranch();
                final NormalSufficientStatistics parent = statistics.getParent();

                if (DEBUG) {
                    System.err.println("B = " + statistics.toVectorizedString());
                }

                DenseMatrix64F Qi = parent.getRawPrecision();
                DenseMatrix64F Si = matrix1; //new DenseMatrix64F(dim, dim); // matrix1;
                CommonOps.scale(differentialScaling, branch.getRawVariance(), Si);

                if (DEBUG) {
                    System.err.println("\tQi = " + NormalSufficientStatistics.toVectorizedString(Qi));
                    System.err.println("\tSi = " + NormalSufficientStatistics.toVectorizedString(Si));
                }

                DenseMatrix64F logDetComponent = matrix0; //new DenseMatrix64F(dim, dim); //matrix0;
                CommonOps.mult(Qi, Si, logDetComponent);

                double grad1 = 0.0;
                for (int row = 0; row < dim; ++row) {
                    grad1 -= 0.5 * logDetComponent.unsafe_get(row, row);
                }

                if (DEBUG) {
                    System.err.println("grad1 = " + grad1);
                }

                DenseMatrix64F quadraticComponent = matrix1; //new DenseMatrix64F(dim, dim); //matrix1;
                CommonOps.mult(logDetComponent, Qi, quadraticComponent);

                if (DEBUG) {
                    System.err.println("\tQi = " + NormalSufficientStatistics.toVectorizedString(Qi));
                    System.err.println("\tQQ = " + NormalSufficientStatistics.toVectorizedString(quadraticComponent));
                }

                NormalSufficientStatistics jointStatistics = computeJointStatistics(child, parent);

                DenseMatrix64F additionalVariance = matrix0; //new DenseMatrix64F(dim, dim);
                CommonOps.mult(quadraticComponent, jointStatistics.getRawVariance(), additionalVariance);

                DenseMatrix64F delta = vector0;
                for (int row = 0; row < dim; ++row) {
                    delta.unsafe_set(row, 0,
                            jointStatistics.getRawMean().unsafe_get(row, 0) - parent.getMean(row)
                    );
                }

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
                    System.err.println("\tjoint mean = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawMean()));
                    System.err.println("\tparent mean = " + NormalSufficientStatistics.toVectorizedString(parent.getRawMean()));
                    System.err.println("\tchild mean = " + NormalSufficientStatistics.toVectorizedString(child.getRawMean()));
                    System.err.println("\tjoint variance = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawVariance()));
                    System.err.println("\tchild variance = " + NormalSufficientStatistics.toVectorizedString(child.getRawVariance()));
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

            public NormalSufficientStatistics computeJointStatistics(NormalSufficientStatistics child,
                                                                      NormalSufficientStatistics parent) {

                DenseMatrix64F totalP = new DenseMatrix64F(dim, dim);
                CommonOps.add(child.getRawPrecision(), parent.getRawPrecision(), totalP);

                DenseMatrix64F totalV = new DenseMatrix64F(dim, dim);
                safeInvert(totalP, totalV, false);

                DenseMatrix64F mean = new DenseMatrix64F(dim, 1);
                safeWeightedAverage(
                        new WrappedVector.Raw(child.getRawMean().getData(), 0, dim),
                        child.getRawPrecision(),
                        new WrappedVector.Raw(parent.getRawMean().getData(), 0, dim),
                        parent.getRawPrecision(),
                        new WrappedVector.Raw(mean.getData(), 0, dim),
                        totalV,
                        dim);

                return new NormalSufficientStatistics(mean, totalP, totalV);
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

    @Override
    public String getReport() {

        double[] savedValues = rateParameter.getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, rateParameter.getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            rateParameter.setParameterValue(i, savedValues[i]);
        }

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

        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("gradient report", new Object() {
            @Override
            public String toString() {
                return "\n" + getReport();
            }
        });

        return columns;
    }
}
