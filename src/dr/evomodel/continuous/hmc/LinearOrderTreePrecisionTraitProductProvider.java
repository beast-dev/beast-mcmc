package dr.evomodel.continuous.hmc;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.WrappedMeanPrecision;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;

import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */
public class LinearOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;
    private static final boolean DEBUG = false;

    public LinearOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                        String traitName) {
        super(treeDataLikelihood, likelihoodDelegate);

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));

        this.delta = new double[tree.getExternalNodeCount()][dimTrait];


//        int threadCount = 0;
//        if (threadCount > 0) {
//            pool = Executors.newFixedThreadPool(threadCount);
//        } else if (threadCount < 0) {
//            // create a cached thread pool which should create one thread per likelihood...
//            pool = Executors.newCachedThreadPool();
//        } else {
//            // don't use a threadpool (i.e., compute serially)
//            pool = null;
//        }

        this.doParallel = false;

    }

    private static final boolean NEW_DATA = false; // Maybe not useful

    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        final double[] result = new double[vector.getDimension()];

        if (!doParallel) { // single-threaded

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


//            try {
//                List<Future<Integer>> results = pool.invokeAll(likelihoodCallers);
//
//                for (Future<Double> result : results) {
//                    double logL = result.get();
//                    logLikelihood += logL;
//                }
//
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            } catch (ExecutionException e) {
//                e.printStackTrace();
//            }


            int threadCount = 4;

            ExecutorService es = Executors.newFixedThreadPool(threadCount);

            for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {

                final int t = taxon;
//
//                    es.execute(new Runnable() {
//                        public void run() {
                            computeProductForOneTaxon(t, allStatistics.get(t), result);
//                        }
//                    });

//                try {
//                    es.submit(new Callable<Integer>() {
//                        public Integer call() {
//                            return computeProductForOneTaxon(t, allStatistics.get(t), result);
//                        }
//                    }).get();
//                } catch (InterruptedException e) {
//                    e.printStackTrace();
//                } catch (ExecutionException e) {
//                    e.printStackTrace();
//                }
            }
        }

        if (DEBUG) {
            debug(result, vector);
        }

        return result;
    }

    private int computeProductForOneTaxon(final int taxon,
                                          final WrappedMeanPrecision statistic,
                                          final double[] result) {

        final ReadableVector mean = statistic.getMean();
        final ReadableMatrix precision = statistic.getPrecision();
        final double scalar = statistic.getPrecisionScalar();
        final int resultOffset = taxon * dimTrait;

        computeDelta(taxon, delta[taxon], dataParameter, mean);
        computePrecisionDeltaProduct(result, resultOffset, precision, delta[taxon], scalar);

        return 1;
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
                sum +=  precision.get(i, j) * delta[j];
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
        return 0.0; // TODO
    }

    @SuppressWarnings("unchecked")
    private TreeTrait<List<WrappedMeanPrecision>> castTreeTrait(TreeTrait trait) {
        return trait;
    }

    private final boolean doParallel;
//    private final ExecutorService pool;
    private final double[][] delta;
}
