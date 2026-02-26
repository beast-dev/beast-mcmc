package dr.inference.model;

public class FiniteMixtureParameter extends Parameter.Proxy implements VariableListener {

    private final Parameter values;
    private final Parameter categories;

    private final int maxCategory;

    private boolean statisticsKnown;

    public FiniteMixtureParameter(String name,
                                  Parameter values,
                                  Parameter categories) {
        super(name, categories.getDimension());

        this.values = values;
        this.categories = categories;
        this.maxCategory = values.getDimension();


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

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        statisticsKnown = false;
    }
}
