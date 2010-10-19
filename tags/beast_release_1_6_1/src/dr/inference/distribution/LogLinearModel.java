package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class LogLinearModel extends GeneralizedLinearModel {

    public LogLinearModel(Parameter dependentParam) {
        super(dependentParam);
    }

    @Override
    public double[] getXBeta() {
        double[] xBeta = super.getXBeta();
        for(int i=0; i<xBeta.length; i++)
            xBeta[i] = Math.exp(xBeta[i]);
        return xBeta;
    }

    protected double calculateLogLikelihood(double[] beta) {
        throw new RuntimeException("Not yet implemented.");
    }

    protected double calculateLogLikelihood() {
        throw new RuntimeException("Not yet implemented.");
    }

    protected boolean confirmIndependentParameters() {
        return false;
    }

    public boolean requiresScale() {
        return false;
    }
}
