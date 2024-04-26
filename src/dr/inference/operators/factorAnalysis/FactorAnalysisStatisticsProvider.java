package dr.inference.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.continuous.HashedMissingArray;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class FactorAnalysisStatisticsProvider implements VariableListener {

    private final FactorAnalysisOperatorAdaptor adaptor;
    private final boolean useInnerProductCache;
    private Map<HashedMissingArray, DenseMatrix64F> innerProductMap = new HashMap<>();

    private final double[][] observedIndicators;
    private boolean needToUpdateCache = true;
    private Parameter[] factorDependentParameters;


    public FactorAnalysisStatisticsProvider(FactorAnalysisOperatorAdaptor adaptor, CacheProvider cacheProvider) {
        this.adaptor = adaptor;
        this.useInnerProductCache = cacheProvider.useCache();

        if (cacheProvider.useCache()) {
            observedIndicators = setupObservedIndicators();

        } else {
            observedIndicators = null;
        }

        this.factorDependentParameters = adaptor.getFactorDependentParameters();
        for (Parameter parameter : factorDependentParameters) {
            parameter.addParameterListener(this);
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

    public void getFactorMeans(double[] buffer) {
        int n = adaptor.getNumberOfTaxa();
        int k = adaptor.getNumberOfFactors();

        for (int factor = 0; factor < k; factor++) {
            double sum = 0;
            for (int taxon = 0; taxon < n; taxon++) {
                sum += adaptor.getFactorValue(factor, taxon);
            }
            buffer[factor] = sum / n;
        }
    }

    public void getLoadingsInnerProduct(double[][] buffer) {
        int k = adaptor.getNumberOfFactors();
        int p = adaptor.getNumberOfTraits();
        for (int k1 = 0; k1 < k; k1++) {
            for (int k2 = k1; k2 < k; k2++) {
                double sum = 0;
                int k1Offset = k1 * p;
                int k2Offset = k2 * p;

                for (int j = 0; j < p; j++) {
                    double l1 = adaptor.getLoadingsValue(j + k1Offset);
                    double l2 = adaptor.getLoadingsValue(j + k2Offset);

                    sum += l1 * l2;
                }

                buffer[k1][k2] = sum;
                if (k1 != k2) {
                    buffer[k2][k1] = sum;
                }
            }
        }

    }

    public void getFactorInnerProduct(int trait, int dim, double[][] buffer) {
        getFactorInnerProduct(trait, dim, buffer, true);
    }

    public void getFactorInnerProduct(int trait, int dim, double[][] buffer, boolean checkMissing) {
        DenseMatrix64F hashedPrecision = null;
        HashedMissingArray observedArray = null;

        if (useInnerProductCache) {
            checkNeedToUpdateInnerProduct();
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
                        if (!checkMissing || adaptor.isNotMissing(trait, k)) {
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

    public void getScaledFactorInnerProduct(int trait, int dim, double[][] buffer) {
        getFactorInnerProduct(trait, dim, buffer);

        double precision = adaptor.getColumnPrecision(trait);

        for (int i = 0; i < dim; i++) {
            buffer[i][i] *= precision;

            for (int j = (i + 1); j < dim; j++) {
                buffer[i][j] *= precision;
                buffer[j][i] = buffer[i][j];
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

    public void getScaledFactorTraitProduct(int trait, int dim, double[] buffer) {
        getFactorTraitProduct(trait, dim, buffer);
        double precision = adaptor.getColumnPrecision(trait);
        for (int i = 0; i < dim; i++) {
            buffer[i] *= precision;
        }
    }

    public void checkNeedToUpdateInnerProduct() {
        if (useInnerProductCache && needToUpdateCache) {
            innerProductMap.clear();
            needToUpdateCache = false;
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        for (Parameter parameter : factorDependentParameters) {
            if (variable == parameter) {
                needToUpdateCache = true;
            }
        }
    }

    public FactorAnalysisOperatorAdaptor getAdaptor() {
        return adaptor;
    } //TODO: maybe make StatisticsProvider extend adaptor?

    public boolean useCache() {
        return useInnerProductCache;
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
