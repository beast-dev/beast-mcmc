package dr.evomodel.branchratemodel;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 */
public class ContinuousEpochBranchRateModel extends RateEpochBranchRateModel {

    /**
     * The constructor. For an N-epoch model, there should be N rate paramters and N-1 transition times.
     *
     * @param timeParameters an array of transition time parameters
     * @param rateParameters an array of rate parameters
     */
    public ContinuousEpochBranchRateModel(Parameter[] timeParameters, Parameter[] rateParameters, Parameter rootHeight) {
        super(timeParameters, rateParameters);
        this.rootHeight = rootHeight;
        addVariable(rootHeight);
        normalizationKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);
        normalizationKnown = false;
    }

    private void normalize() {
        normalization = 0.0;
        double startTime = 0.0;
        double endTime = rootHeight.getParameterValue(0);
        int j = 0;
        while( j < timeParameters.length && endTime > timeParameters[j].getParameterValue(0)) {
            final double nextTime = timeParameters[j].getParameterValue(0);
            normalization += (nextTime - startTime) * rateParameters[j].getParameterValue(0);
            startTime = nextTime;
            j++;
        }
        normalization += (endTime - startTime) * rateParameters[j].getParameterValue(0);
    }

    protected void storeState() {
        savedNormalization = normalization;
    }

    protected void restoreState() {
        normalization = savedNormalization;
    }

    protected double normalizeRate(double rate) {
        if (!normalizationKnown)
            normalize();
        return rate / normalization;
    }

    private Parameter rootHeight;
    private double normalization;
    private double savedNormalization;
    private boolean normalizationKnown = false;

}
