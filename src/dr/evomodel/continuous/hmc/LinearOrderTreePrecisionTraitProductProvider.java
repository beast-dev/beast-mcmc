package dr.evomodel.continuous.hmc;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.NewTipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class LinearOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    private final TreeTrait<List<NormalSufficientStatistics>> fullConditionalDensity;
    private static final boolean DEBUG = false;

    public LinearOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                        String traitName) {
        super(treeDataLikelihood, likelihoodDelegate);

        String fcdName = NewTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addNewFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
    }

    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        double[] result = new double[vector.getDimension()];

        int offset = 0;
        for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {
            List<NormalSufficientStatistics> statistics = fullConditionalDensity.getTrait(
                    tree, tree.getExternalNode(taxon));

            assert (statistics.size() == 1);

            NormalSufficientStatistics statistic = statistics.get(0);
            for (int i = 0; i < dimTrait; ++i) {
                double sum = 0.0;
                for (int j = 0; j < dimTrait; ++j) {
                    sum += statistic.getPrecision(i, j) *
                            (dataParameter.getParameterValue(taxon * dimTrait + j) - statistic.getMean(j));
                }
                result[offset] = sum;
                ++offset;
            }
        }

        if (DEBUG) {

            double[] result2 = expensiveProduct(vector);

            System.err.println("via FCD: " + new WrappedVector.Raw(result));
            System.err.println("direct : " + new WrappedVector.Raw(result2));
            System.err.println();

        }

        return result;
    }

    @Override
    public double[] getMassVector() {
        return null; // TODO
    }

    @Override
    public double getTimeScale() {
        return 0.0; // TODO
    }

    @SuppressWarnings("unchecked")
    private TreeTrait<List<NormalSufficientStatistics>> castTreeTrait(TreeTrait trait) {
        return trait;
    }
}
