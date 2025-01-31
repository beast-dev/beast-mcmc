package dr.inference.model;

import java.util.List;

/**
 * This is parameter class that is intended to allow sampling values on the range
 * [0,) inclusive of 0. It represents a specific case of a produce parameter.

 */
public class IndicatorProductParameter extends Parameter.Abstract implements VariableListener {
    Parameter indicatorParameter;
    Parameter continuousParameter;
    
    private Bounds bounds = null;

    public IndicatorProductParameter(Parameter indicatorParameter, Parameter continuousParameter) {

        this.indicatorParameter = indicatorParameter;
        this.continuousParameter = continuousParameter;

        indicatorParameter.addVariableListener(this);
        Parameter.CONNECTED_PARAMETER_SET.add(indicatorParameter);

        continuousParameter.addVariableListener(this);
        Parameter.CONNECTED_PARAMETER_SET.add(continuousParameter);
    }
    public int getDimension() {
        return indicatorParameter.getDimension();
    }
    @Override
    public boolean isImmutable() {
        return true;
    }

    protected void storeValues() {
        indicatorParameter.storeParameterValues();
        continuousParameter.storeParameterValues();
    }
    protected void restoreValues() {
        indicatorParameter.restoreParameterValues();
        continuousParameter.restoreParameterValues();
    }
    protected void acceptValues() {
        indicatorParameter.acceptParameterValues();
        continuousParameter.acceptParameterValues();
    }

    public double getParameterValue(int dim) {
        double value = indicatorParameter.getParameterValue(dim) *continuousParameter.getParameterValue(dim);
        return value;
    }

    public void setParameterValue(int dim, double value) {
        if(value==0){
            indicatorParameter.setParameterValue(dim, 0);
        }else{
            indicatorParameter.setParameterValue(dim, 1);
            continuousParameter.setParameterValue(dim, value);
        }
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }
    @Override
    protected void adoptValues(Parameter source) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'adoptValues'");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("IndicatorProduct");
           
            sb.append(".").append(indicatorParameter.getId());
            sb.append(".").append(continuousParameter.getId());

            setId(sb.toString());
        }
        return getId();
    }
    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        if (bounds == null) {
            return indicatorParameter.getBounds(); // TODO
        } else {
            return bounds;
        }
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index,type);
    }



}
