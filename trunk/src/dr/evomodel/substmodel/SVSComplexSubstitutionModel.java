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

    protected void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        if (parameter == infinitesimalRates && indicators.getParameterValue(index) == 0)
            return;  // ignore, does not affect likelihood
        super.handleParameterChangedEvent(parameter, index, type);
    }

    public Parameter getIndicators() {
        return indicators;
    }

     public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0) { // Also check that graph is connected
            getTransitionProbabilities(1.0,testProbabilities);
            for(double value : testProbabilities) {
                if (value <= 0 || Double.isInfinite(value) || Double.isNaN(value)) {
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
