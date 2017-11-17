package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.NormalDistribution;

public class IndependentNormalDistributionModel extends AbstractModelLikelihood {
    Parameter mean;
    Parameter variance;
    Parameter precision;
    Parameter data;
    boolean usePrecision;

    public IndependentNormalDistributionModel(String id, Parameter mean, Parameter variance, Parameter precision, Parameter data){
        super(id);
        addVariable(mean);
        this.mean = mean;
        if(precision != null){
            usePrecision = true;
            addVariable(precision);
        }
        else{
            usePrecision = false;
            addVariable(variance);
        }
        this.precision = precision;
        this.variance = variance;
        this.data = data;
        addVariable(data);
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double sum = 0;
        for (int i = 0; i < data.getDimension(); i++) {
            double sd;
            if(usePrecision){
                sd = Math.sqrt(1 / precision.getParameterValue(i));
            }
            else{
                sd = Math.sqrt(variance.getParameterValue(i));
            }
            sum += NormalDistribution.logPdf(data.getParameterValue(i),
                    mean.getParameterValue(i),
                    sd);
        }
        return sum;
    }

    @Override
    public void makeDirty() {

    }

    public Parameter getPrecision() {
        return precision;
    }

    public Parameter getVariance() {
        return variance;
    }

    public Parameter getMean() {
        return mean;
    }
}
