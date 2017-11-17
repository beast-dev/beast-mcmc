package dr.inference.distribution;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.InverseGammaDistribution;

public class IndependentInverseGammaDistributionModel extends AbstractModelLikelihood {
        Parameter shape;
        Parameter scale;
        Parameter data;

    public IndependentInverseGammaDistributionModel(String id, Parameter shape, Parameter scale, Parameter data){
        super(id);
        this.shape = shape;
        addVariable(shape);
        this.scale = scale;
        addVariable(scale);
        this.data = data;
        addVariable(data);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

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
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double sum = 0;
        for (int i = 0; i < data.getDimension(); i++) {
            sum += InverseGammaDistribution.logPdf(data.getParameterValue(i),
                    shape.getParameterValue(i), scale.getParameterValue(i), 0.0);
        }
        return sum;
    }

    public Parameter getShape() {
        return shape;
    }

    public Parameter getScale() {
        return scale;
    }

    public Parameter getData() {
        return data;
    }

    @Override
    public void makeDirty() {

    }
}
