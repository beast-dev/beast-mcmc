package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

public interface OUGradientStrategy {

    DenseMatrix64F getGradientVarianceWrtVariance(OUDiffusionModelDelegate delegate,
                                                  NodeRef node,
                                                  ContinuousDiffusionIntegrator cdi,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                  DenseMatrix64F gradient);

    DenseMatrix64F getGradientVarianceWrtAttenuation(OUDiffusionModelDelegate delegate,
                                                     NodeRef node,
                                                     ContinuousDiffusionIntegrator cdi,
                                                     BranchSufficientStatistics statistics,
                                                     DenseMatrix64F gradient);

    DenseMatrix64F getGradientDisplacementWrtAttenuation(OUDiffusionModelDelegate delegate,
                                                         NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F gradient);

    DenseMatrix64F getGradientDisplacementWrtDrift(OUDiffusionModelDelegate delegate,
                                                   NodeRef node,
                                                   ContinuousDiffusionIntegrator cdi,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   DenseMatrix64F gradient);

    double[] getGradientDisplacementWrtRoot(OUDiffusionModelDelegate delegate,
                                            NodeRef node,
                                            ContinuousDiffusionIntegrator cdi,
                                            ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                            DenseMatrix64F gradient);
}
