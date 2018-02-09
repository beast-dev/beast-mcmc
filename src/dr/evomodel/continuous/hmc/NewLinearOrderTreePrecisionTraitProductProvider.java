package dr.evomodel.continuous.hmc;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.WrappedMeanPrecision;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class NewLinearOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;
    private static final boolean DEBUG = false;

    public NewLinearOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                           ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                           String traitName) {
        super(treeDataLikelihood, likelihoodDelegate);

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));

        this.delta = new double[dimTrait];
    }

    private static final boolean NEW_COMPUTE = true;
    private static final boolean NEW_DATA = true;

    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        double[] result = new double[vector.getDimension()];

        List<WrappedMeanPrecision> statistics;
        int statisticIndex = 0;

        if (NEW_DATA) {
            statistics = fullConditionalDensity.getTrait(tree, null);
            assert (statistics.size() == tree.getExternalNodeCount());
        }

        int offset = 0;
        for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) { // TODO In parallel

            if (!NEW_DATA) {
                statistics = fullConditionalDensity.getTrait(tree, tree.getExternalNode(taxon)); // TODO Get once
                assert (statistics.size() == 1);
            }
            
            final WrappedMeanPrecision statistic = statistics.get(statisticIndex);
            final ReadableVector mean = statistic.getMean();
            final ReadableMatrix precision = statistic.getPrecision();
            final double scalar = statistic.getPrecisionScalar();

            if (NEW_COMPUTE) {

                computeDelta(taxon, delta, dataParameter, mean);
                computePrecisionDeltaProduct(result, offset, precision, delta, scalar);
                offset += dimTrait;

            } else {

                for (int i = 0; i < dimTrait; ++i) {
                    delta[i] = dataParameter.getParameterValue(taxon * dimTrait + i) - mean.get(i);
                }

                for (int i = 0; i < dimTrait; ++i) {
                    double sum = 0.0;
                    for (int j = 0; j < dimTrait; ++j) {
                        sum += scalar * precision.get(i, j) * delta[j];
                    }
                    result[offset] = sum;
                    ++offset;
                }
                
            }

            if (NEW_DATA) {
                ++statisticIndex;
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

    private static void computeDelta(final int taxon,
                                     final double[] delta,
                                     final Parameter data,
                                     final ReadableVector mean) {
        for (int i = 0, dim = delta.length; i < dim; ++i) {
            delta[i] = data.getParameterValue(taxon * dim + i) - mean.get(i);
        }
    }

    private static void computePrecisionDeltaProduct(final double[] result,
                                                     int offset,
                                                     final ReadableMatrix precision,
                                                     final double[] delta,
                                                     final double scalar) {
        final int dim = delta.length;

        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dim; ++j) {
                sum += scalar * precision.get(i, j) * delta[j];
            }
            result[offset] = sum;
            ++offset;
        }
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
    private TreeTrait<List<WrappedMeanPrecision>> castTreeTrait(TreeTrait trait) {
        return trait;
    }

    private final double[] delta;
}
