package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.evolution.datatype.DataType;

/**
 * @author Marc Suchard
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements BayesianStochasticSearchVariableSelection {

    public SVSGeneralSubstitutionModel(DataType dataType, FrequencyModel freqModel,
                                       Parameter ratesParameter, Parameter indicatorsParameter) {
        super(dataType, freqModel, ratesParameter, -1);
        this.indicatorsParameter  = indicatorsParameter;
        addVariable(indicatorsParameter);
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++) {
            rates[i] = ratesParameter.getParameterValue(i) * indicatorsParameter.getParameterValue(i);
        }
    }

    public Parameter getIndicators() {
        return indicatorsParameter;
    }

    public boolean validState() {
        return true;
    }

    private final Parameter indicatorsParameter;

}
