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

    private final boolean useInnerProductCache = true; //TODO: in parser
    private Map<HashedMissingArray, DenseMatrix64F> factorInnerProductMap = new HashMap<>();

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

            double[][] fullAnswer = answer;
            int rowDimension = newRowDimension;
            if (useInnerProductCache) {
                fullAnswer = new double[nFactors][nFactors];
                rowDimension = nFactors;
            }

            int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

            for (int i = 0; i < rowDimension; i++) {
                for (int j = i; j < rowDimension; j++) {
                    double sum = 0;
                    for (int k = 0; k < p; k++) {
                        if (adaptor.isNotMissing(trait, k)) {
                            sum += adaptor.getFactorValue(i, k) * adaptor.getFactorValue(j, k);
                        }
                    }

                    fullAnswer[i][j] = sum;
                    if (i != j) {
                        fullAnswer[j][i] = sum;
                    }
                }
            }

            if (useInnerProductCache) {

                for (int i = 0; i < newRowDimension; i++) {
                    System.arraycopy(fullAnswer[i], 0, answer[i], 0, newRowDimension);
                }

                factorInnerProductMap.put(observedArray, new DenseMatrix64F(fullAnswer));
            }
        } else {
            for (int i = 0; i < newRowDimension; i++) {
                System.arraycopy(hashedInnerProduct.getData(), i * nFactors,
                        answer[i], 0, newRowDimension);
            }

        }
    }


    public void getFactorTraitProductForTrait(int trait, int newRowDimension, double[] answer) {

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


    }


}
