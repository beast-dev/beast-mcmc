package dr.evomodel.continuous.hmc;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.WrappedMeanPrecision;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;

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

    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;

    private static final boolean DEBUG = false;
    private static final boolean NEW_DATA = false; // Maybe not useful
    private static final boolean SMART_POOL = true;

    public LinearOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                        String traitName,
                                                        int threadCount,
                                                        double roughTimeGuess,
                                                        int eigenvalueReplicates) {
        super(treeDataLikelihood, likelihoodDelegate);

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));

        this.delta = new double[tree.getExternalNodeCount()][dimTrait];
		
        this.roughTimeGuess = roughTimeGuess;
        this.eigenvalueReplicates = eigenvalueReplicates;
        
        setupParallelServices(tree.getExternalNodeCount(), threadCount);
    }


    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        final double[] result = new double[vector.getDimension()];

        if (pool == null) { // single-threaded

            final List<WrappedMeanPrecision> allStatistics;
            if (NEW_DATA) {
                allStatistics = fullConditionalDensity.getTrait(tree, null);
                assert (allStatistics.size() == tree.getExternalNodeCount());
            }

            for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {

                final WrappedMeanPrecision statistic;
                if (NEW_DATA) {
                    statistic = allStatistics.get(taxon);
                } else {
                    List<WrappedMeanPrecision> statistics = fullConditionalDensity.getTrait(tree, tree.getExternalNode(taxon));
                    assert (statistics.size() == 1);
                    statistic = statistics.get(0);
                }

                computeProductForOneTaxon(taxon, statistic, result);
            }

        } else {

            final List<WrappedMeanPrecision> allStatistics = fullConditionalDensity.getTrait(tree, null);
            assert (allStatistics.size() == tree.getExternalNodeCount());

            List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

            if (SMART_POOL) {

                for (final TaskIndices indices : taskIndices) {
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
                pool.invokeAll(calls);
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
                                          final WrappedMeanPrecision statistic,
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

        if (roughTimeGuess > 0.0) { // TODO Super bad, some delegate for re-use with other Providers
            return roughTimeGuess;
        }
        return getRoughLowerBoundforTravelTime();
    }
    
    private double getMaxEigenvalueAsTravelTime() {

        // TODO Lots of bad magic numbers
        return maxEigenvalueByPowerMethod(likelihoodDelegate.getTreeVariance(), 50, 0.01)
                *  maxEigenvalueByPowerMethod(likelihoodDelegate.getTraitVariance(), 50, 0.01);
    }

    private static double maxEigenvalueByPowerMethod(double[][] matrix, int numIterations, double err) {

        double[] y0 = new double[matrix.length];
        ReadableVector diff;
        double maxEigenvalue = 10.0; // TODO Bad magic number
 
        for (int i = 0; i < matrix.length; ++i) {
            y0[i] = MathUtils.nextDouble();
        }
        WrappedVector y = new WrappedVector.Raw(y0);

        final ReadableMatrix mat = new WrappedMatrix.ArrayOfArray(matrix);

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
        return maxEigenvalue;
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
    private TreeTrait<List<WrappedMeanPrecision>> castTreeTrait(TreeTrait trait) {
        return trait;
    }

    private void setupParallelServices(int taxonCount, int threadCount) {
        if (threadCount > 0) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            pool = null;
        }

        taskIndices = (pool != null) ? setupTasks(taxonCount, threadCount) : null;
    }

    private List<TaskIndices> setupTasks(int taxonCount, int threadCount) {
        List<TaskIndices> tasks = new ArrayList<TaskIndices>(threadCount);

        int length = taxonCount / threadCount;
        if (taxonCount % threadCount != 0) ++length;

        int start = 0;

        for (int task = 0; task < threadCount && start < taxonCount; ++task) {
            tasks.add(new TaskIndices(start, Math.min(start + length, taxonCount)));
            start += length;
        }

        return tasks;
    }

    private class TaskIndices {
        int start;
        int stop;

        TaskIndices(int start, int stop) {
            this.start = start;
            this.stop = stop;
        }

        public String toString() {
            return start + " " + stop;
        }
    }

    private ExecutorService pool = null;
    private List<TaskIndices> taskIndices = null;

    private final double[][] delta;
    private final double roughTimeGuess;
    private final int eigenvalueReplicates;
}
