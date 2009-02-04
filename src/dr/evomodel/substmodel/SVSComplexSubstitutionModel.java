package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class SVSComplexSubstitutionModel extends ComplexSubstitutionModel implements BayesianStochasticSearchVariableSelection {

    public SVSComplexSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel, Parameter parameter, Parameter indicators) {
        super(name, dataType, rootFreqModel, parameter);
        this.indicators = indicators;
        addParameter(indicators);
        testProbabilities = new double[stateCount*stateCount];
    }

    protected double[] getRates() {
        double[] rates = infinitesimalRates.getParameterValues();
        double[] indies = indicators.getParameterValues();
        final int dim = rates.length;
        for (int i = 0; i < dim; i++)
            rates[i] *= indies[i];
        return rates;
    }

    public Parameter getIndicators() {
        return indicators;
    }

     public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0) { // Also check that graph is connected
            getTransitionProbabilities(1.0,testProbabilities);
            for(int i=0; i<testProbabilities.length; i++) {
                if (testProbabilities[i] == 0) {
                    logL = Double.NEGATIVE_INFINITY;
                    break;
                }
            }
        }
        return logL;
    }

    public boolean validState() {
        return true;
    }

    private Parameter indicators;
    private double[] testProbabilities;

}
