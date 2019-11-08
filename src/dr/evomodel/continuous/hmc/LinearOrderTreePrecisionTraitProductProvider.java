package dr.evomodel.continuous.hmc;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

import static dr.math.matrixAlgebra.ReadableMatrix.Utils.product;
import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;
import static dr.math.matrixAlgebra.ReadableVector.Utils.norm;

/**
 * @author Marc A. Suchard
 */
public class LinearOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    private final TreeTrait<List<WrappedNormalSufficientStatistics>> fullConditionalDensity;

    private static final boolean DEBUG = false;
    private static final boolean NEW_DATA = false; // Maybe not useful
    private static final boolean SMART_POOL = true;

    public LinearOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                        String traitName,
                                                        int threadCount,
                                                        double roughTimeGuess,
                                                        double optimalTravelTimeScalar,
                                                        int eigenvalueReplicates) {
        super(treeDataLikelihood, likelihoodDelegate);

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));

        this.delta = new double[tree.getExternalNodeCount()][dimTrait];
		
        this.roughTimeGuess = roughTimeGuess;
        this.optimalTravelTimeScalar = optimalTravelTimeScalar;
        this.eigenvalueReplicates = eigenvalueReplicates;

        this.taxonTaskPool = new TaxonTaskPool(tree.getExternalNodeCount(), threadCount);
    }
    
    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        final double[] result = new double[vector.getDimension()];

        if (taxonTaskPool.getPool() == null) { // single-threaded

            final List<WrappedNormalSufficientStatistics> allStatistics;
            if (NEW_DATA) {
                allStatistics = fullConditionalDensity.getTrait(tree, null);
                assert (allStatistics.size() == tree.getExternalNodeCount());
            }

            for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {

                final WrappedNormalSufficientStatistics statistic;
                if (NEW_DATA) {
                    statistic = allStatistics.get(taxon);
                } else {
                    List<WrappedNormalSufficientStatistics> statistics = fullConditionalDensity.getTrait(tree, tree.getExternalNode(taxon));
                    assert (statistics.size() == 1);
                    statistic = statistics.get(0);
                }

                computeProductForOneTaxon(taxon, statistic, result);
            }

        } else {

            final List<WrappedNormalSufficientStatistics> allStatistics = fullConditionalDensity.getTrait(tree, null);
            assert (allStatistics.size() == tree.getExternalNodeCount());

            List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

            if (SMART_POOL) {

                for (final TaxonTaskPool.TaxonTaskIndices indices : taxonTaskPool.getIndices()) {
                    calls.add(Executors.callable(
                            new Runnable() {
                                @Override
                                public void run() {
                                    for (int taxon = indices.start; taxon < indices.stop; ++taxon) {
                                        computeProductForOneTaxon(taxon, allStatistics.get(taxon), result);
                                    }
                                }
                            }
                    ));
                }

            } else {

                for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {

                    final int t = taxon;
                    calls.add(Executors.callable(
                            new Runnable() {
                                public void run() {
                                    computeProductForOneTaxon(t, allStatistics.get(t), result);
                                }
                            }
                    ));
                }
            }

            try {
                taxonTaskPool.getPool().invokeAll(calls);
            } catch (InterruptedException exception) {
                exception.printStackTrace();
            }
        }

        if (DEBUG) {
            debug(result, vector);
        }

        return result;
    }

    private void computeProductForOneTaxon(final int taxon,
                                           final WrappedNormalSufficientStatistics statistic,
                                           final double[] result) {

        final ReadableVector mean = statistic.getMean();
        final ReadableMatrix precision = statistic.getPrecision();
        final double scalar = statistic.getPrecisionScalar();
        final int resultOffset = taxon * dimTrait;

        computeDelta(taxon, delta[taxon], dataParameter, mean);
        computePrecisionDeltaProduct(result, resultOffset, precision, delta[taxon], scalar);
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
                sum += precision.get(i, j) * delta[j];
            }
            result[offset] = sum * scalar;
            ++offset;
        }
    }

    @Override
    public double[] getMassVector() {
        return null; // TODO
    }

    @Override
    public double getTimeScale() {

        if (roughTimeGuess > 0.0) {
            return roughTimeGuess;
        } else {
            return getMaxEigenvalueAsTravelTime();
        }
    }

    @Override
    public double getTimeScaleEigen() {
        return maxEigenvalueByPowerMethod(likelihoodDelegate.getTraitVariance(), 50, 0.01, false); //TODO: magic numbers
    }

    private double getMaxEigenvalueAsTravelTime() {

        // TODO Lots of bad magic numbers
        double treeCovEigenValue = maxEigenvalueByPowerMethod(likelihoodDelegate.getTreeVariance(), 50, 0.01, false);
        double traitCovEigenValue = maxEigenvalueByPowerMethod(likelihoodDelegate.getTraitVariance(), 50, 0.01, false);
        return optimalTravelTimeScalar * Math.sqrt(treeCovEigenValue * traitCovEigenValue);
    }

    private static double maxEigenvalueByPowerMethod(double[][] matrix, int numIterations, double err, boolean inverseflag) {

        double[][] matrixForUse;

        if (inverseflag) {
            matrixForUse = (new SymmetricMatrix(matrix)).inverse().toComponents();
        } else {
            matrixForUse = (new SymmetricMatrix(matrix)).toComponents();
        }

        double[] y0 = new double[matrixForUse.length];
        ReadableVector diff;
        double maxEigenvalue = 10.0; // TODO Bad magic number
 
        for (int i = 0; i < matrixForUse.length; ++i) {
            y0[i] = MathUtils.nextDouble();
        }
        WrappedVector y = new WrappedVector.Raw(y0);

        final ReadableMatrix mat = new WrappedMatrix.ArrayOfArray(matrixForUse);

        for (int i = 0; i < numIterations; ++i) {

            ReadableVector v = new ReadableVector.Scale(1 / norm(y), y);
            y = product(mat, v);
            maxEigenvalue = innerProduct(v, y);
            diff = new ReadableVector.Sum(y,
                    new ReadableVector.Scale(-maxEigenvalue, v));

            if (ReadableVector.Utils.norm(diff) < err) {
                break;
            }
        }

        if (inverseflag) {
            return 1.0 / maxEigenvalue;
        } else {
            return maxEigenvalue;
        }

    }

    private double getRoughLowerBoundforTravelTime() {

        ReadableVector savedDataParameter = new WrappedVector.Raw(dataParameter.getParameterValues());

        double precisionMinEigenvalueLowerBound = 0.0;
        for (int i = 0; i < eigenvalueReplicates; ++i) {

            ReadableVector x = drawUniformSphere(dataParameter.getDimension());
            ReadableVector.Utils.setParameter(x, dataParameter);

            ReadableVector Phi_x = new WrappedVector.Raw(getProduct(dataParameter));

            precisionMinEigenvalueLowerBound += ReadableVector.Utils.innerProduct(x, Phi_x);

        }
        precisionMinEigenvalueLowerBound /= eigenvalueReplicates; // TODO Could compute average on sqrt(1/bound) scale

        ReadableVector.Utils.setParameter(savedDataParameter, dataParameter);

        return Math.sqrt(1 / precisionMinEigenvalueLowerBound);
    }

    private static WrappedVector drawUniformSphere(final int len) {

        double[] x = new double[len];
        double normSquare = 0.0;

        for (int i = 0; i < len; i++) {
            x[i] = MathUtils.nextGaussian();
            normSquare += x[i] * x[i];
        }

        double norm = Math.sqrt(normSquare);

        for (int i = 0; i < len; i++) {
            x[i] = x[i] / norm;
        }

        return new WrappedVector.Raw(x);
    }

    @SuppressWarnings("unchecked")
    static TreeTrait<List<WrappedNormalSufficientStatistics>> castTreeTrait(TreeTrait trait) {
        return trait;
    }

    private final TaxonTaskPool taxonTaskPool;

    private final double[][] delta;
    private final double roughTimeGuess;
    private final int eigenvalueReplicates;
    private final double optimalTravelTimeScalar;
}
