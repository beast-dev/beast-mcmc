package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class ComplexSubstitutionModel extends GeneralSubstitutionModel {

    public ComplexSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, Parameter parameter) {
        super(name, dataType, freqModel, parameter, -1);
    }

    protected EigenSystem getDefaultEigenSystem(int stateCount) {
        return new ColtEigenSystem();
    }

    protected int getRateCount(int stateCount) {
        return (stateCount - 1) * stateCount;
    }

    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++)
            rates[i] = ratesParameter.getParameterValue(i);
    }

    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
