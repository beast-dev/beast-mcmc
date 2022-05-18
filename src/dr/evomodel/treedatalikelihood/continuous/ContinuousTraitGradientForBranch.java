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
import dr.evomodel.treedatalikelihood.hmc.MultivariateChainRule;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public interface ContinuousTraitGradientForBranch {

    double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node);

//    double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node, boolean getGradientQ, boolean getGradientN);

    int getParameterIndexFromNode(NodeRef node);

    int getDimension();

    abstract class Default implements ContinuousTraitGradientForBranch {

        private final DenseMatrix64F matrixGradientQInv;
        private final DenseMatrix64F matrixGradientN;

        final DenseMatrix64F matrixDelta;
        DenseMatrix64F matrixQ;
        DenseMatrix64F matrixW;
        DenseMatrix64F matrixV;

        final int dim;
        final Tree tree;

        public Default(int dim, Tree tree) {
            this.dim = dim;
            this.tree = tree;

            matrixGradientQInv = new DenseMatrix64F(dim, dim);
            matrixGradientN = new DenseMatrix64F(dim, 1);

            matrixDelta = new DenseMatrix64F(dim, 1);
            matrixQ = new DenseMatrix64F(dim, dim);
            matrixW = new DenseMatrix64F(dim, dim);
            matrixV = new DenseMatrix64F(dim, dim);
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return node.getNumber();
        }

        @Override
        public double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node) {
            return getGradientForBranch(statistics, node, true, true);
        }

        double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node,
                                      boolean getGradientQ, boolean getGradientN) {

            getSufficientStatistics(statistics, node);

            DenseMatrix64F Qi = matrixQ;
            DenseMatrix64F Wi = matrixW;
            DenseMatrix64F Vi = matrixV;
            DenseMatrix64F delta = matrixDelta;

            DenseMatrix64F gradQInv = matrixGradientQInv;
            DenseMatrix64F gradN = matrixGradientN;
            if (getGradientQ) getGradientQInvForBranch(Qi, Wi, Vi, delta, gradQInv);
            if (getGradientN) getGradientNForBranch(Qi, delta, gradN);

            if (tree.isRoot(node)) {
                return chainRuleRoot(statistics, node, gradQInv, gradN);
            } else {
                return chainRule(statistics, node, gradQInv, gradN);
            }
        }

        void getSufficientStatistics(BranchSufficientStatistics statistics, NodeRef node) {
            // Joint Statistics
            final NormalSufficientStatistics below = statistics.getBelow();
            final NormalSufficientStatistics above = statistics.getAbove();
            NormalSufficientStatistics jointStatistics =
                    BranchRateGradient.ContinuousTraitGradientForBranch.Default.computeJointStatistics(
                            below, above, dim
                    );

            matrixQ = above.getRawPrecision();
            matrixW = above.getRawVariance();
            matrixV = jointStatistics.getRawVariance();

            if (DEBUG) {
                System.err.println("B = " + statistics.toVectorizedString());
                System.err.println("\tjoint mean = " + NormalSufficientStatistics.toVectorizedString(jointStatistics.getRawMean()));
                System.err.println("\tabove mean = " + NormalSufficientStatistics.toVectorizedString(above.getRawMean()));
                System.err.println("\tbelow mean = " + NormalSufficientStatistics.toVectorizedString(below.getRawMean()));
                System.err.println("\tjoint variance Vi = " + NormalSufficientStatistics.toVectorizedString(matrixV));
                System.err.println("\tbelow variance = " + NormalSufficientStatistics.toVectorizedString(below.getRawVariance()));
                System.err.println("\tabove variance Wi = " + NormalSufficientStatistics.toVectorizedString(matrixW));
                System.err.println("\tabove precision Qi = " + NormalSufficientStatistics.toVectorizedString(matrixQ));
            }

            // Delta
            for (int row = 0; row < dim; ++row) {
                matrixDelta.unsafe_set(row, 0,
                        jointStatistics.getMean(row) - above.getMean(row)
                );
            }

            if (DEBUG) {
                System.err.println("\tDelta = " + NormalSufficientStatistics.toVectorizedString(matrixDelta));
            }
        }

        abstract double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                                    DenseMatrix64F gradQInv, DenseMatrix64F gradN);

        abstract double[] chainRuleRoot(BranchSufficientStatistics statistics, NodeRef node,
                                        DenseMatrix64F gradQInv, DenseMatrix64F gradN);

        static void getGradientQInvForBranch(DenseMatrix64F Qi, DenseMatrix64F Wi,
                                             DenseMatrix64F Vi, DenseMatrix64F delta,
                                             DenseMatrix64F grad) {

            CommonOps.scale(0.5, Wi, grad);

            CommonOps.multAddTransB(-0.5, delta, delta, grad);

            CommonOps.addEquals(grad, -0.5, Vi);

            if (DEBUG) {
                System.err.println("\tgradientQi = " + NormalSufficientStatistics.toVectorizedString(grad));
            }

            MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(Qi);
            ruleI.chainGradient(grad);

            if (DEBUG) {
                System.err.println("\tgradientQiInv = " + NormalSufficientStatistics.toVectorizedString(grad));
            }

        }

        private void getGradientNForBranch(DenseMatrix64F Qi, DenseMatrix64F delta, DenseMatrix64F grad) {

            CommonOps.multTransA(Qi, delta, grad);

            if (DEBUG) {
                System.err.println("\tgradientNi = " + NormalSufficientStatistics.toVectorizedString(grad));
            }

        }

        static final boolean DEBUG = false;
    }

    class RateGradient extends Default {

        DiffusionProcessDelegate diffusionProcessDelegate;

        private final DenseMatrix64F matrixJacobianQInv;
        private final DenseMatrix64F matrixJacobianN;
        private final DenseMatrix64F matrix0;

        private final ArbitraryBranchRates branchRateModel;

        public RateGradient(int dim, Tree tree, BranchRateModel brm,
                            DiffusionProcessDelegate diffusionProcessDelegate) {
            super(dim, tree);

            this.branchRateModel = (brm instanceof ArbitraryBranchRates) ? (ArbitraryBranchRates) brm : null;

            matrixJacobianQInv = new DenseMatrix64F(dim, dim);
            matrixJacobianN = new DenseMatrix64F(dim, 1);
            matrix0 = new DenseMatrix64F(dim, dim);

            this.diffusionProcessDelegate = diffusionProcessDelegate;
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
                                  DenseMatrix64F gradQInv, DenseMatrix64F gradN) {

            final double rate = branchRateModel.getBranchRate(tree, node);
            final double differential = branchRateModel.getBranchRateDifferential(tree, node);
            final double scaling = differential / rate;

            // Q_i w.r.t. rate
            DenseMatrix64F gradMatQInv = matrixJacobianQInv;
            CommonOps.scale(scaling, statistics.getBranch().getRawVariance(), gradMatQInv);

            double[] gradient = new double[1];
            for (int i = 0; i < gradMatQInv.getNumElements(); i++) {
                gradient[0] += gradMatQInv.get(i) * gradQInv.get(i);
            }

            // n_i w.r.t. rate
            DenseMatrix64F gradMatN = matrixJacobianN;
            diffusionProcessDelegate.updateGradientDisplacementWrtRate(gradient, scaling,
                    statistics.getBranch().getRawDisplacement(),
                    gradMatN, gradN);

            return gradient;

        }

        @Override
        public double[] chainRuleRoot(BranchSufficientStatistics statistics, NodeRef node,
                                      DenseMatrix64F gradQInv, DenseMatrix64F gradN) {

            return new double[1];

        }
    }

    class ContinuousProcessParameterGradient extends Default {

        ContinuousDataLikelihoodDelegate likelihoodDelegate;
        ContinuousDiffusionIntegrator cdi;
        DiffusionProcessDelegate diffusionProcessDelegate;

        final List<DerivationParameter> derivationParameter;

        public ContinuousProcessParameterGradient(int dim, Tree tree,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                  List<DerivationParameter> derivationParameter) {
            super(dim, tree);

            this.likelihoodDelegate = likelihoodDelegate;
            this.cdi = likelihoodDelegate.getIntegrator();
            this.diffusionProcessDelegate = likelihoodDelegate.getDiffusionProcessDelegate();
            this.derivationParameter = derivationParameter;

        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return 0;
        }

        @Override
        public int getDimension() {
            int paramDim = 0;
            for (DerivationParameter param : derivationParameter) {
                paramDim += param.getDimension(dim);
            }
            return paramDim;
        }

        @Override
        public double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                                  DenseMatrix64F gradQInv, DenseMatrix64F gradN) {

//            cdi.getVariancePreOrderDerivative(statistics, gradQ);
            removeMissing(gradQInv, statistics.getMissing());

            double[] gradient = new double[getDimension()];
            int offset = 0;
            for (DerivationParameter param : derivationParameter) {
                int paramDim = param.getDimension(dim);
                System.arraycopy(
                        param.chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN),
                        0, gradient, offset, paramDim);
                offset += paramDim;
            }
            return gradient;
        }

        @Override
        public double[] chainRuleRoot(BranchSufficientStatistics statistics, NodeRef node,
                                      DenseMatrix64F gradQInv, DenseMatrix64F gradN) {

            double[] gradient = new double[getDimension()];
            int offset = 0;
            for (DerivationParameter param : derivationParameter) {
                int paramDim = param.getDimension(dim);
                System.arraycopy(
                        param.chainRuleRoot(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN),
                        0, gradient, offset, paramDim);
                offset += paramDim;
            }
            return gradient;
        }

        private static void removeMissing(DenseMatrix64F M, int[] missing) {
            for (int m : missing) {
                for (int j = 0; j < M.getNumCols(); j++) {
                    M.unsafe_set(m, j, 0.0);
                    M.unsafe_set(j, m, 0.0);
                }
            }
        }

        List<DerivationParameter> getDerivationParameter() {
            return derivationParameter;
        }

        public enum DerivationParameter {
            WRT_VARIANCE {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    DenseMatrix64F gradient = diffusionProcessDelegate.getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradQInv);

                    if (DEBUG) {
                        System.err.println("gradQ = " + NormalSufficientStatistics.toVectorizedString(gradient));
                    }

                    return gradient.getData();

                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    return chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);

                }

                @Override
                public int getDimension(int dim) {
                    return dim * dim;
                }
            },
            WRT_CONSTANT_DRIFT {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    DenseMatrix64F gradient = ((AbstractDriftDiffusionModelDelegate) diffusionProcessDelegate).getGradientDisplacementWrtDrift(node, cdi, likelihoodDelegate, gradN);

                    if (DEBUG) {
                        System.err.println("gradQ = " + NormalSufficientStatistics.toVectorizedString(gradient));
                    }

                    return gradient.getData();
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    return new double[gradN.getNumRows()];
                }

                @Override
                public int getDimension(int dim) {
                    return dim;
                }
            },
            WRT_ROOT_MEAN {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    return diffusionProcessDelegate.getGradientDisplacementWrtRoot(node, cdi, likelihoodDelegate, gradN);
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    return chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
                }

                @Override
                public int getDimension(int dim) {
                    return dim;
                }
            },
            WRT_CONSTANT_DRIFT_AND_ROOT_MEAN {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    double[] drift = WRT_CONSTANT_DRIFT.chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
                    double[] root = WRT_ROOT_MEAN.chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
                    for (int i = 0; i < root.length; i++) {
                        root[i] += drift[i];
                    }
                    return root;
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    return WRT_ROOT_MEAN.chainRuleRoot(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
                }

                @Override
                public int getDimension(int dim) {
                    return dim;
                }
            },
            WRT_DIAGONAL_SELECTION_STRENGTH {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    DenseMatrix64F gradQInvDiag = ((OUDiffusionModelDelegate) diffusionProcessDelegate).getGradientVarianceWrtAttenuation(node, cdi, statistics, gradQInv);

                    if (DEBUG) {
                        System.err.println("gradQ = " + NormalSufficientStatistics.toVectorizedString(gradQInv));
                    }

                    DenseMatrix64F gradNDiag = ((OUDiffusionModelDelegate) diffusionProcessDelegate).getGradientDisplacementWrtAttenuation(node, cdi, statistics, gradN);

                    CommonOps.addEquals(gradQInvDiag, gradNDiag);

                    return gradQInvDiag.getData();
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              DenseMatrix64F gradQInv, DenseMatrix64F gradN) {
                    return new double[likelihoodDelegate.getTraitDim()];
                }

                @Override
                public int getDimension(int dim) {
                    return dim;
                }
            },

            WRT_BRANCH_SPECIFIC_DRIFT {

                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {

                    double[] drift = WRT_CONSTANT_DRIFT.chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
//                    double[] root = WRT_ROOT_MEAN.chainRule(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
//                    for (int i = 0; i < root.length; i++) {
//                        root[i] += drift[i];
//                    }
                    // return root
                    return drift;
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    return WRT_ROOT_MEAN.chainRuleRoot(cdi, diffusionProcessDelegate, likelihoodDelegate, statistics, node, gradQInv, gradN);
                }

                @Override
                public int getDimension(int dim) {
                    return dim;
                }
            },
            WRT_SAMPLING_VARIANCE {
                @Override
                public double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                          DiffusionProcessDelegate diffusionProcessDelegate,
                                          ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                          BranchSufficientStatistics statistics, NodeRef node,
                                          final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    return gradQInv.getData();
                }

                @Override
                public double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                              DiffusionProcessDelegate diffusionProcessDelegate,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              BranchSufficientStatistics statistics, NodeRef node,
                                              final DenseMatrix64F gradQInv, final DenseMatrix64F gradN) {
                    throw new RuntimeException("Should never be called.");
                }

                @Override
                public int getDimension(int dim) {
                    return dim * dim;
                }
            };

            abstract double[] chainRule(ContinuousDiffusionIntegrator cdi,
                                        DiffusionProcessDelegate diffusionProcessDelegate,
                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                        BranchSufficientStatistics statistics, NodeRef node,
                                        final DenseMatrix64F gradQInv, final DenseMatrix64F gradN);

            abstract double[] chainRuleRoot(ContinuousDiffusionIntegrator cdi,
                                            DiffusionProcessDelegate diffusionProcessDelegate,
                                            ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                            BranchSufficientStatistics statistics, NodeRef node,
                                            final DenseMatrix64F gradQInv, final DenseMatrix64F gradN);

            public abstract int getDimension(int dim);
        }
    }

    class SamplingVarianceGradient extends ContinuousProcessParameterGradient {

        ModelExtensionProvider.NormalExtensionProvider dataModel;


        public SamplingVarianceGradient(int dim,
                                        Tree tree,
                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                        ModelExtensionProvider.NormalExtensionProvider dataModel) {
            super(dim, tree,
                    likelihoodDelegate,
                    Arrays.asList(DerivationParameter.WRT_SAMPLING_VARIANCE));

            this.dataModel = dataModel;
        }

        @Override
        public double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node) {
            if (!tree.isExternal(node)) {
                return new double[getDimension()]; // 0 if not a tip
            }

            double[] gradient = getGradientForBranch(statistics, node, true, false);
            dataModel.chainRuleWrtVariance(gradient, node);
            return gradient;
        }

        void getSufficientStatistics(BranchSufficientStatistics statistics, NodeRef node) {
            final NormalSufficientStatistics below = statistics.getBelow();
            final NormalSufficientStatistics above = statistics.getAbove();

            // Sampling Parameters: Statistic data model
            DenseMatrix64F samplingVariance = dataModel.getExtensionVariance(node);

            // One more pre-order step
            // TODO This is just one more pre-order step. Should maybe be moved elsewhere ?
            // TODO Below only works for one data point (no repetition)
            // above data
            DenseMatrix64F aboveDataVariance = new DenseMatrix64F(dim, dim);
            DenseMatrix64F aboveDataPrecision = new DenseMatrix64F(dim, dim);
            CommonOps.add(above.getRawVariance(), samplingVariance, aboveDataVariance);
            safeInvert2(aboveDataVariance, aboveDataPrecision, false);
            final NormalSufficientStatistics aboveData = new NormalSufficientStatistics(
                    above.getRawMeanCopy(),
                    aboveDataPrecision,
                    aboveDataVariance);
            // Below data
            int[] missings = statistics.getMissing();
            DenseMatrix64F precisionData = new DenseMatrix64F(dim, dim);
            for (int i = 0; i < dim; i++) {
                precisionData.unsafe_set(i, i, Double.POSITIVE_INFINITY);
            }
            for (int m : missings) {
                precisionData.unsafe_set(m, m, 0.0);
            }
            final NormalSufficientStatistics belowData = new NormalSufficientStatistics(
                    below.getRawMeanCopy(),
                    precisionData);

            // Joint data
            NormalSufficientStatistics jointStatisticsData =
                    BranchRateGradient.ContinuousTraitGradientForBranch.Default.computeJointStatistics(
                            belowData, aboveData, dim
                    );

            // Gradient
            matrixQ = aboveData.getRawPrecision();
            matrixW = aboveData.getRawVariance();
            matrixV = jointStatisticsData.getRawVariance();

            // Delta
            for (int row = 0; row < dim; ++row) {
                matrixDelta.unsafe_set(row, 0,
                        jointStatisticsData.getMean(row) - aboveData.getMean(row)
                );
            }
        }

        @Override
        public int getParameterIndexFromNode(NodeRef node) {
            return 0;
        }

    }

}
