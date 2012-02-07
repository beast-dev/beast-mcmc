package dr.inference.model;

/**
 * @author Marc Suchard
 */
public class DuplicatedParameter extends Parameter.Abstract implements VariableListener {

    public DuplicatedParameter(Parameter parameter) {
        this.parameter = parameter;
        parameter.addVariableListener(this);
        copies = 1;
        originalBounds = parameter.getBounds();
        bounds = originalBounds;
    }

    public void addDuplicationParameter(Parameter dupParameter) {
        this.dupParameter = dupParameter;
        dupParameter.addParameterListener(this);
        updateDuplication();
    }

    private void updateDuplication() {
        copies = (int) dupParameter.getParameterValue(0);
        final int originalLength = originalBounds.getBoundsDimension();
        double[] lowers = new double[getDimension()];
        double[] uppers = new double[getDimension()];
        for(int i=0; i<originalLength; i++) {
            lowers[i] = originalBounds.getLowerLimit(i);
            uppers[i] = originalBounds.getUpperLimit(i);
        }
        for(int i=1; i<copies; i++) {
            System.arraycopy(lowers,0,lowers,i*originalLength,originalLength);
            System.arraycopy(uppers,0,uppers,i*originalLength,originalLength);
        }
        bounds = new DefaultBounds(uppers,lowers);
    }

    public int getDimension() {
        return parameter.getDimension() * copies;
    }

    protected void storeValues() {
        parameter.storeParameterValues();
    }

    protected void restoreValues() {
        parameter.restoreParameterValues();
    }

    protected void acceptValues() {
        parameter.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        parameter.adoptParameterValues(source);
    }

    public double getParameterValue(int dim) {
        return parameter.getParameterValue(dim % parameter.getDimension());
    }

    public void setParameterValue(int dim, double value) {
        parameter.setParameterValue(dim % parameter.getDimension(), value);
        fireParameterChangedEvent(dim, Parameter.ChangeType.VALUE_CHANGED);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(dim % parameter.getDimension(), value);
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        parameter.setParameterValueNotifyChangedAll(dim % parameter.getDimension(), value);
    }

    public String getParameterName() {
        if (getId() == null)
            return "duplicated" + parameter.getParameterName();
        return getId();
    }

    public void addBounds(Bounds bounds) {
        throw new RuntimeException("Not yet implemented.");
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if (variable == dupParameter) {
            updateDuplication();
        }
        // Values have changed, so notify listeners
        for (int i = 0; i < copies; ++i) {
            // fire once for each duplication
            fireParameterChangedEvent(index + i * parameter.getDimension(), Parameter.ChangeType.VALUE_CHANGED);
        }
    }

    private final Parameter parameter;
    private Parameter dupParameter;
    private int copies;
    private Bounds<Double> bounds;
    private final Bounds<Double> originalBounds;

}
