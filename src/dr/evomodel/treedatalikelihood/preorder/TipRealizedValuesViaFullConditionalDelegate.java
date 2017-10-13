package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 */
public class TipRealizedValuesViaFullConditionalDelegate extends AbstractValuesViaFullConditionalDelegate {

    public TipRealizedValuesViaFullConditionalDelegate(String name, Tree tree,
                                                       MultivariateDiffusionModel diffusionModel,
                                                       ContinuousTraitPartialsProvider dataModel,
                                                       ConjugateRootTraitPrior rootPrior,
                                                       ContinuousRateTransformation rateTransformation,
//                                                       BranchRateModel rateModel,
                                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public static String getName(String name) {
        return "tipSample." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

    @Override
    protected void computeValueWithNoMissing(final double[] mean, final int meanOffset,
                                             final double[] output, final int outputOffset,
                                             final int dim) {
        System.arraycopy(mean, meanOffset, output, outputOffset, dim);
    }

    @Override
    protected void computeValueWithMissing(final WrappedVector mean,
                                           final double[][] cholesky,
                                           final WrappedVector output,
                                           final double[] buffer) {
        MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                mean, // input mean
                cholesky, 1.0, // input variance
                output, // output sample
                buffer);
    }
}
