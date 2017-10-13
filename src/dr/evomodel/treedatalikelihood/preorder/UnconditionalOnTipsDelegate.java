package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.math.distributions.MultivariateNormalDistribution;

/**
 * @author Marc A. Suchard
 */
public class UnconditionalOnTipsDelegate extends AbstractRealizedContinuousTraitDelegate {

    public UnconditionalOnTipsDelegate(String name,
                                       MutableTreeModel tree,
                                       MultivariateDiffusionModel diffusionModel,
                                       ContinuousTraitDataModel dataModel,
                                       ConjugateRootTraitPrior rootPrior,
                                       ContinuousRateTransformation rateTransformation,
                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        this.rootPrior = rootPrior;
    }

    @Override
    protected void simulateRoot(final int nodeIndex) {

        final double[] rootMean = rootPrior.getMean();
        final double sqrtScale = Math.sqrt(1.0 / rootPrior.getPseudoObservations());

        int offsetSample = dimNode * nodeIndex;
        for (int trait = 0; trait < numTraits; ++trait) {
            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    rootMean, 0, // input mean
                    cholesky, sqrtScale,
                    sample, offsetSample,
                    tmpEpsilon
            );

            offsetSample += dimTrait;
        }
    }

    protected void simulateNode(final BranchNodeOperation operation,
                                final double branchNormalization) {
        final int nodeIndex = operation.getNodeNumber();
        int offsetSample = dimNode * nodeIndex;
        int offsetParent = dimNode * operation.getParentNumber();

        final double branchLength = operation.getBranchLength() * branchNormalization;

        if (branchLength == 0.0) {
            System.arraycopy(sample, offsetParent, sample, offsetSample, dimTrait * numTraits);
        } else {

            final double sqrtScale = Math.sqrt(branchLength);
            for (int trait = 0; trait < numTraits; ++trait) {
                MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                        sample, offsetParent,
                        cholesky, sqrtScale,
                        sample, offsetSample,
                        tmpEpsilon
                );

                offsetParent += dimTrait;
                offsetSample += dimTrait;
            }
        }
    }

    @Override
    protected void simulateNode(final int parentNumber,
                                         final int nodeNumber,
                                         final int nodeMatrix,
                                         final int siblingNumber,
                                         final int siblingMatrix) {
        throw new RuntimeException("Not yet implemented -- see above");
    }

    private final ConjugateRootTraitPrior rootPrior;
}
