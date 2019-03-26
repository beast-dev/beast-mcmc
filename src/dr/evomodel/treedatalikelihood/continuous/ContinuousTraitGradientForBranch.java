/*
 * ContinuousTraitGradientForBranch.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert;
import static dr.math.matrixAlgebra.missingData.MissingOps.safeWeightedAverage;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public interface ContinuousTraitGradientForBranch {

    double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node);

    double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node, boolean getGradientQ, boolean getGradientN);

    int getParameterIndexFromNode(NodeRef node);

    int getDimension();

    abstract class Default implements ContinuousTraitGradientForBranch {

        private final DenseMatrix64F matrixGradientQ;
        private final DenseMatrix64F matrixGradientN;
        private final DenseMatrix64F vector0;

        final int dim;
        final Tree tree;

        public Default(int dim, Tree tree) {
            this.dim = dim;
            this.tree = tree;

            matrixGradientQ = new DenseMatrix64F(dim, dim);
            matrixGradientN = new DenseMatrix64F(dim, 1);
            vector0 = new DenseMatrix64F(dim, 1);
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return node.getNumber();
        }

        @Override
        public double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node) {
            return getGradientForBranch(statistics, node, true, true);
        }

        @Override
        public double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node,
                                             boolean getGradientQ, boolean getGradientN) {
            // Joint Statistics
            final NormalSufficientStatistics child = statistics.getChild();
            final NormalSufficientStatistics parent = statistics.getParent();
            NormalSufficientStatistics jointStatistics =
                    BranchRateGradient.ContinuousTraitGradientForBranch.Default.computeJointStatistics(
                            child, parent, dim
                    );
//                    computeJointStatistics(child, parent);

            DenseMatrix64F Qi = parent.getRawPrecision();
            DenseMatrix64F Wi = parent.getRawVariance();
            DenseMatrix64F Vi = jointStatistics.getRawVariance();

            if (DEBUG) {
                System.err.println("B = " + statistics.toVectorizedString());
                System.err.println("\tjoint mean = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawMean()));
                System.err.println("\tparent mean = " + NormalSufficientStatistics.toVectorizedString(parent.getRawMean()));
                System.err.println("\tchild mean = " + NormalSufficientStatistics.toVectorizedString(child.getRawMean()));
                System.err.println("\tjoint variance Vi = " + NormalSufficientStatistics.toVectorizedString(Vi));
                System.err.println("\tchild variance = " + NormalSufficientStatistics.toVectorizedString(child.getRawVariance()));
                System.err.println("\tparent variance Wi = " + NormalSufficientStatistics.toVectorizedString(Wi));
                System.err.println("\tparent precision Qi = " + NormalSufficientStatistics.toVectorizedString(Qi));
            }

            // Delta
            DenseMatrix64F delta = vector0;
            for (int row = 0; row < dim; ++row) {
                delta.unsafe_set(row, 0,
                        jointStatistics.getRawMean().unsafe_get(row, 0) - parent.getMean(row)
                );
            }

            if (DEBUG) {
                System.err.println("\tDelta = " + NormalSufficientStatistics.toVectorizedString(delta));
            }

            DenseMatrix64F gradQ = matrixGradientQ;
            DenseMatrix64F gradN = matrixGradientN;
            if (getGradientQ) getGradientQForBranch(Wi, Vi, delta, gradQ);
            if (getGradientN) getGradientNForBranch(Qi, delta, gradN);

            if (tree.isRoot(node)) {
                return chainRuleRoot(node, gradQ, gradN);
            } else {
                return chainRule(statistics, node, gradQ, gradN);
            }
        }

        abstract double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                                    DenseMatrix64F gradQ, DenseMatrix64F gradN);

        abstract double[] chainRuleRoot(NodeRef node, DenseMatrix64F gradQ, DenseMatrix64F gradN);

        private void getGradientQForBranch(DenseMatrix64F Wi, DenseMatrix64F Vi, DenseMatrix64F delta,
                                           DenseMatrix64F grad) {

            CommonOps.scale(0.5, Wi, grad);

            CommonOps.multAddTransB(-0.5, delta, delta, grad);

            CommonOps.addEquals(grad, -0.5, Vi);

            if (DEBUG) {
                System.err.println("\tgradientQi = " + NormalSufficientStatistics.toVectorizedString(grad));
            }

        }

        private void getGradientNForBranch(DenseMatrix64F Qi, DenseMatrix64F delta, DenseMatrix64F grad) {

            CommonOps.multTransA(Qi, delta, grad);

            if (DEBUG) {
                System.err.println("\tgradientNi = " + NormalSufficientStatistics.toVectorizedString(grad));
            }

        }

        private NormalSufficientStatistics computeJointStatistics(NormalSufficientStatistics child,
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

        static final boolean DEBUG = false;
    }

    class RateGradient extends Default {

        private final DenseMatrix64F matrixJacobianQ;
        private final DenseMatrix64F matrixJacobianN;
        private final DenseMatrix64F matrix0;

        private final ArbitraryBranchRates branchRateModel;

        public RateGradient(int dim, Tree tree, BranchRateModel brm) {
            super(dim, tree);

            this.branchRateModel = (brm instanceof ArbitraryBranchRates) ? (ArbitraryBranchRates) brm : null;

            matrixJacobianQ = new DenseMatrix64F(dim, dim);
            matrixJacobianN = new DenseMatrix64F(dim, 1);
            matrix0 = new DenseMatrix64F(dim, dim);
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            if (tree.isRoot(node)) return 0;
            return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
        }

        @Override
        public int getDimension() {
            return 1;
        }

        @Override
        public double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                                  DenseMatrix64F gradQ, DenseMatrix64F gradN) {

            final double rate = branchRateModel.getBranchRate(tree, node);
            final double differential = branchRateModel.getBranchRateDifferential(tree, node);
            final double scaling = differential / rate;

            // Q_i w.r.t. rate
            DenseMatrix64F gradMatQ = matrixJacobianQ;
            CommonOps.scale(scaling, statistics.getBranch().getRawVariance(), gradMatQ);

            DenseMatrix64F Qi = statistics.getParent().getRawPrecision();
            DenseMatrix64F temp = matrix0;
            CommonOps.mult(Qi, gradMatQ, temp);
            CommonOps.mult(-1.0, temp, Qi, gradMatQ);

            double[] gradVecQ = gradMatQ.getData();

            double[] gradient = new double[1];
            for (int i = 0; i < gradVecQ.length; i++) {
                gradient[0] += gradVecQ[i] * gradQ.get(i);
            }

            // n_i w.r.t. rate
            // TODO: Fix delegate to (possibly) un-link drift from arbitrary rate
            DenseMatrix64F gradMatN = matrixJacobianN;
            CommonOps.scale(scaling, statistics.getBranch().getRawDisplacement(), gradMatN);
            for (int i = 0; i < gradMatN.numRows; i++) {
                gradient[0] += gradMatN.get(i) * gradN.get(i);
            }

            return gradient;

        }

        @Override
        public double[] chainRuleRoot(NodeRef node, DenseMatrix64F gradQ, DenseMatrix64F gradN) {

            return new double[1];

        }
    }

    class ContinuousProcessParameterGradient extends Default {
        final DenseMatrix64F matrixJacobianQ;
        final DenseMatrix64F matrixJacobianN;
        final DenseMatrix64F matrix0;

        ContinuousDataLikelihoodDelegate likelihoodDelegate;
        ContinuousDiffusionIntegrator cdi;
        DiffusionProcessDelegate diffusionProcessDelegate;

        final DerivationParameter derivationParameter;

        public ContinuousProcessParameterGradient(int dim, Tree tree,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                  DerivationParameter derivationParameter) {
            super(dim, tree);

            this.likelihoodDelegate = likelihoodDelegate;
            this.cdi = likelihoodDelegate.getIntegrator();
            this.diffusionProcessDelegate = likelihoodDelegate.getDiffusionProcessDelegate();
            this.derivationParameter = derivationParameter;

            matrixJacobianQ = new DenseMatrix64F(dim, dim);
            matrixJacobianN = new DenseMatrix64F(dim, 1);
            matrix0 = new DenseMatrix64F(dim, dim);
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return 0;
        }

        @Override
        public int getDimension() {
            return dim * dim;
        }

        public enum DerivationParameter {
            WRT_PRECISION {
                @Override
                public void preOrderGradientPrecision(ContinuousDiffusionIntegrator cdi,
                                                      BranchSufficientStatistics statistics, DenseMatrix64F gradientQ) {
                    cdi.getPrecisionPreOrderDerivative(statistics, gradientQ);
                }
            },
            WRT_DRIFT {
                @Override
                public void preOrderGradientPrecision(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics, DenseMatrix64F gradientQ) {
                    throw new RuntimeException("not yet implemented");
                }
            },
            WRT_SELECTION_STRENGTH {
                @Override
                public void preOrderGradientPrecision(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics, DenseMatrix64F gradientQ) {
                    throw new RuntimeException("not yet implemented");
                }
            };

            abstract void preOrderGradientPrecision(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics, DenseMatrix64F gradientQ);
        }

        @Override
        public double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                                  DenseMatrix64F gradQ, DenseMatrix64F gradN) {

            derivationParameter.preOrderGradientPrecision(cdi, statistics, gradQ);

            return chainRuleRoot(node, gradQ, gradN);

        }

        @Override
        public double[] chainRuleRoot(NodeRef node, DenseMatrix64F gradQ, DenseMatrix64F gradN) {

            diffusionProcessDelegate.getGradientPrecision(getScalarNode(node), gradQ);

            if (DEBUG) {
                System.err.println("gradQ = " + NormalSufficientStatistics.toVectorizedString(gradQ));
            }

            return gradQ.getData();

        }

//        private int getPreOrderPartialIndex(NodeRef node) {
//            return likelihoodDelegate.getPartialBufferCount() + node.getNumber();
//        }

        private int getActiveMatrixIndex(NodeRef node) {
            return likelihoodDelegate.getActiveMatrixIndex(node.getNumber());
        }

        private double getScalarNode(NodeRef node) {
            if (tree.isRoot(node)) {
                return likelihoodDelegate.getRootProcessDelegate().getPseudoObservations();
            } else {
                return cdi.getInverseBranchLength(getActiveMatrixIndex(node));
            }
        }
    }

}
