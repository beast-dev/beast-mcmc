package dr.inference.model;

public class FiniteMixtureParameter extends Parameter.Proxy implements VariableListener {

    private final Parameter values;
    private final Parameter categories;

    private final int maxCategory;

    public FiniteMixtureParameter(String name,
                                  Parameter values,
                                  Parameter categories) {
        super(name, categories.getDimension());

        this.values = values;
        this.categories = categories;
        this.maxCategory = values.getDimension();

        Parameter.CONNECTED_PARAMETER_SET.add(values);
        Parameter.CONNECTED_PARAMETER_SET.add(categories);

        values.addParameterListener(this);
        categories.addParameterListener(this);
    }

    public int[] getOccupancy() {
        int[] occupancy = new int[maxCategory];
        for (int i = 0; i < categories.getDimension(); ++i) {
            int category = (int) categories.getParameterValue(i);
            ++occupancy[category];
        }
        return occupancy;
    }

    public Parameter getValuesParameter() { return values; }

    public Parameter getCategoriesParameter() { return categories; }

    @Override
    public double getParameterValue(int dim) {
        int category = (int) categories.getParameterValue(dim);
        if (category >= maxCategory) {
            throw new IllegalArgumentException("Invalid category");
        }
        return values.getParameterValue(category);
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public Bounds<Double> getBounds() {
        throw new RuntimeException("Do not operate directly on finite-mixture parameter '" + getId() + "'");
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent();
    }

    @Override
    protected void storeValues() {
        values.storeParameterValues();
        categories.storeParameterValues();
    }

    @Override
    protected void restoreValues() {
        values.restoreParameterValues();
        categories.restoreParameterValues();
    }

    @Override
    protected void acceptValues() {
        values.acceptParameterValues();
        categories.acceptParameterValues();
    }
}
