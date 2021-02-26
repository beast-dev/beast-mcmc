package dr.inference.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.continuous.HashedMissingArray;
import org.ejml.data.DenseMatrix64F;

import java.util.HashMap;
import java.util.Map;

public class LoadingsStatisticsProvider {

    private final FactorAnalysisOperatorAdaptor adaptor;
    private final int nTaxa;
    private final int nTraits;
    private final int nFactors;

    private final boolean useInnerProductCache = false; //TODO: in parser
    private Map<HashedMissingArray, DenseMatrix64F> factorInnerProductMap = new HashMap<>();
    private Map<HashedMissingArray, double[]> factorTraitProductMap = new HashMap<>();

    private static int numThreads = 1; //TODO: in parser
    private static boolean multiThreaded = false; //TODO in parser


    private boolean needToDrawFactors = true; //TODO: turn on and off with listeners (maybe should be in adaptor)
    private final double[][] observedIndicators;


    LoadingsStatisticsProvider(FactorAnalysisOperatorAdaptor adaptor) {
        this.adaptor = adaptor;
        this.nTaxa = adaptor.getNumberOfTaxa();
        this.nTraits = adaptor.getNumberOfTraits();
        this.nFactors = adaptor.getNumberOfFactors();

        if (useInnerProductCache) {
            if (multiThreaded && numThreads > 1) { //TODO: I don't think I need any of this...
                throw new IllegalArgumentException("Cannot currently parallelize cached precisions");
            }
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


    public void getFactorInnerProductForTrait(int trait, int newRowDimension, double[][] answer) {

        DenseMatrix64F hashedInnerProduct = null;
        HashedMissingArray observedArray = null;

        if (useInnerProductCache) {
            double[] observed = observedIndicators[trait];
            observedArray = new HashedMissingArray(observed);
            hashedInnerProduct = factorInnerProductMap.get(observedArray);
        }

        if (!useInnerProductCache || hashedInnerProduct == null) {

            int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

            for (int i = 0; i < newRowDimension; i++) {
                for (int j = i; j < newRowDimension; j++) {
                    double sum = 0;
                    for (int k = 0; k < p; k++) {
                        if (adaptor.isNotMissing(trait, k)) {
                            sum += adaptor.getFactorValue(i, k) * adaptor.getFactorValue(j, k);
                        }
                    }

                    answer[i][j] = sum;
                    if (i != j) {
                        answer[j][i] = sum;
                    }
                }
            }

            if (useInnerProductCache) {
                factorInnerProductMap.put(observedArray, new DenseMatrix64F(answer));
            }
        } else {
            for (int i = 0; i < newRowDimension; i++) {
                System.arraycopy(hashedInnerProduct.getData(), i * newRowDimension,
                        answer[i], 0, newRowDimension);
            }

        }
    }



    public void getFactorTraitProductForTrait(int trait, int newRowDimension, double[] answer) {
        double[] hashedProduct = null;
        HashedMissingArray observedArray = null;

        if (useInnerProductCache) { //TODO: remove code duplication with getFactorInnerProductForTrait?
            double[] observed = observedIndicators[trait];
            observedArray = new HashedMissingArray(observed);
            hashedProduct = factorTraitProductMap.get(observedArray);
        }

        if (!useInnerProductCache || hashedProduct == null) {

            int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

            for (int i = 0; i < newRowDimension; i++) {
                double sum = 0;

                for (int k = 0; k < p; k++) {
                    if (adaptor.isNotMissing(trait, k)) {
                        sum += adaptor.getFactorValue(i, k) /*Left.getParameterValue(i, k)*/
                                * adaptor.getDataValue(trait, k); //data.getParameterValue(dataColumn, k);
                    }
                }


                answer[i] = sum;
            }

            if (useInnerProductCache) {
                double[] answerCopy = new double[newRowDimension];
                System.arraycopy(answer, 0, answerCopy, 0, newRowDimension);
                factorTraitProductMap.put(observedArray, answerCopy);
            }
        } else {
            System.arraycopy(hashedProduct, 0, answer, 0, newRowDimension);
        }
    }



}
