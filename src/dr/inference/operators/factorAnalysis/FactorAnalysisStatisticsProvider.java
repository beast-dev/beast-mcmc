package dr.inference.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.continuous.HashedMissingArray;
import org.ejml.data.DenseMatrix64F;

import java.util.HashMap;
import java.util.Map;

public class FactorAnalysisStatisticsProvider {

    private final FactorAnalysisOperatorAdaptor adaptor;
    private final boolean useInnerProductCache;
    private Map<HashedMissingArray, DenseMatrix64F> innerProductMap = new HashMap<>();

    private final double[][] observedIndicators;


    FactorAnalysisStatisticsProvider(FactorAnalysisOperatorAdaptor adaptor, CacheProvider cacheProvider) {
        this.adaptor = adaptor;
        this.useInnerProductCache = cacheProvider.useCache();

        if (cacheProvider.useCache()) {
            observedIndicators = setupObservedIndicators();

        } else {
            observedIndicators = null;
        }


    }

    private double[][] setupObservedIndicators() {
        double[][] obsInds = new double[adaptor.getNumberOfTraits()][adaptor.getNumberOfTaxa()];

        for (int trait = 0; trait < adaptor.getNumberOfTraits(); trait++) {
            for (int taxon = 0; taxon < adaptor.getNumberOfTaxa(); taxon++) {

                if (adaptor.isNotMissing(trait, taxon)) {
                    obsInds[trait][taxon] = 1.0;
                }
            }
        }


        return obsInds;
    }

    public void getFactorInnerProduct(int trait, int dim, double[][] buffer) {
        DenseMatrix64F hashedPrecision = null;
        HashedMissingArray observedArray = null;

        if (useInnerProductCache) {
            double[] observed = observedIndicators[trait];
            observedArray = new HashedMissingArray(observed);
            hashedPrecision = innerProductMap.get(observedArray);
        }

        if (!useInnerProductCache || hashedPrecision == null) {

            int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    double sum = 0;
                    for (int k = 0; k < p; k++) {
                        if (adaptor.isNotMissing(trait, k)) {
                            sum += adaptor.getFactorValue(i, k) * adaptor.getFactorValue(j, k);
                        }
                    }

                    buffer[i][j] = sum;
                    if (i != j) {
                        buffer[j][i] = sum;
                    }
                }
            }

            if (useInnerProductCache) {
                innerProductMap.put(observedArray, new DenseMatrix64F(buffer));
            }
        } else {
            for (int i = 0; i < dim; i++) {
                System.arraycopy(hashedPrecision.getData(), i * dim,
                        buffer[i], 0, dim);
            }

        }

    }

    public void getFactorTraitProduct(int trait, int dim, double[] buffer) {
        final int p = adaptor.getNumberOfTaxa();

        for (int i = 0; i < dim; i++) {
            double sum = 0;

            for (int k = 0; k < p; k++) {
                if (adaptor.isNotMissing(trait, k)) {
                    sum += adaptor.getFactorValue(i, k) /*Left.getParameterValue(i, k)*/
                            * adaptor.getDataValue(trait, k); //data.getParameterValue(dataColumn, k);
                }
            }

            buffer[i] = sum;
        }
    }

    public void clearInnerProductMap() { //TODO: base on listeners and remove function
        if (useInnerProductCache) {
            innerProductMap.clear();
        }
    }


    public enum CacheProvider {
        USE_CACHE {
            @Override
            boolean useCache() {
                return true;
            }

        },
        NO_CACHE {
            @Override
            boolean useCache() {
                return false;
            }
        };

        abstract boolean useCache();

    }
}
