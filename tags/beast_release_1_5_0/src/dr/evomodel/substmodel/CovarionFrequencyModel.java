package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 */
public class CovarionFrequencyModel extends FrequencyModel {

    public CovarionFrequencyModel(DataType dataType, Parameter frequencyParameter, Parameter hiddenFrequencies) {
        super(dataType, frequencyParameter);

        this.hiddenFrequencies = hiddenFrequencies;
        addVariable(hiddenFrequencies);
    }

    public double[] getFrequencies() {

        int k = 0;

        int numStates = frequencyParameter.getDimension();
        int numHiddenStates = hiddenFrequencies.getDimension();

        double[] freqs = new double[numStates * numHiddenStates];
        for (int i = 0; i < numHiddenStates; i++) {
            for (int j = 0; j < numStates; j++) {
                freqs[k] = frequencyParameter.getParameterValue(j) *
                        hiddenFrequencies.getParameterValue(i);
                k += 1;
            }
        }

        return freqs;
    }

    public void setFrequency(int i, double value) {
        throw new UnsupportedOperationException();
    }

    public double getFrequency(int i) {

        int numStates = frequencyParameter.getDimension();

        return frequencyParameter.getParameterValue(i % numStates) *
                hiddenFrequencies.getParameterValue(i / numStates);
    }

    public int getFrequencyCount() {
        return frequencyParameter.getDimension() * hiddenFrequencies.getDimension();
    }

    Parameter hiddenFrequencies;
}
