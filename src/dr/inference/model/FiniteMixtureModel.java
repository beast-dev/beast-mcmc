package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.distributions.Distribution;
import dr.math.distributions.GeometricDistribution;
import dr.math.distributions.PoissonDistribution;

import java.util.Arrays;

public class FiniteMixtureModel extends AbstractModelLikelihood implements Loggable {

    private final FiniteMixtureParameter mixtureParameter;
    private final Distribution prior;
    private final Constraint constraint;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    public FiniteMixtureModel(String name,
                              FiniteMixtureParameter mixtureParameter,
                              Distribution prior,
                              Constraint constraint) {
        super(name);

        this.mixtureParameter = mixtureParameter;
        this.prior = prior;
        this.constraint = constraint;

        addVariable(mixtureParameter);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Should not get called");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
    }

    @Override
    protected void acceptState() { }

    @Override
    public Model getModel() { return this; }

    private double calculateLikelihood() {

        if (constraint == Constraint.ORDERED) {
            Parameter values = mixtureParameter.getValuesParameter();
            for (int i = 1; i < values.getDimension(); ++i) {
                if (values.getParameterValue(i) - values.getParameterValue(i - 1) <= 0.0) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        int[] occupancy = mixtureParameter.getOccupancy();
        int nonZeroCount = countNonZero(occupancy);
        if (prior instanceof PoissonDistribution) {
            return prior.logPdf(nonZeroCount - 1);
        } else if (prior instanceof GeometricDistribution) {
            return prior.logPdf(nonZeroCount);
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    private int countNonZero(int[] occupancy) {
        int count = 0;
        for (int o : occupancy) {
            if (o > 0) {
                ++count;
            }
        }
        return count;
    }

    @Override
    public double getLogLikelihood() {

        if (!likelihoodKnown) {
            logLikelihood = calculateLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() { likelihoodKnown = false; }

    @Override
    public LogColumn[] getColumns() {

        LogColumn[] superColumns = super.getColumns();

        final Parameter values = mixtureParameter.getValuesParameter();
        final Parameter categories = mixtureParameter.getCategoriesParameter();

        LogColumn[] columns = new LogColumn[superColumns.length + values.getDimension() + categories.getDimension()];

        int offset = 0;
        for (LogColumn column : superColumns) {
            columns[offset++] = column;
        }

        final Integer[] map = new Integer[mixtureParameter.getValuesParameter().getDimension()];

        for (int i = 0; i < values.getDimension(); ++i) {
            columns[offset++] = new SortedValueColumn(
                    values.getParameterName() + ".sorted",
                    values, map, i);
        }

        for (int i = 0; i < categories.getDimension(); ++i) {
            columns[offset++] = new SortedIndexColumn(
                    categories.getParameterName() + ".sorted",
                    categories, map, i);
        }

        return columns;
    }

    private static abstract class SortedColumn extends NumberColumn {

        protected final Parameter parameter;
        protected final Integer[] map;
        protected final int dim;

        public SortedColumn(String label, Parameter parameter, Integer[] map, int dim) {
            super(label);

            this.parameter = parameter;
            this.map = map;
            this.dim = dim;
        }
    }

    private static class SortedIndexColumn extends SortedColumn {

        public SortedIndexColumn(String label, Parameter parameter, Integer[] map, int dim) {
            super(label, parameter, map, dim);
        }

        @Override
        public double getDoubleValue() {
            return map[(int) parameter.getParameterValue(dim)];
        }
    }

    private static class SortedValueColumn extends SortedColumn {

        public SortedValueColumn(String label, Parameter parameter, Integer[] map, int dim) {
            super(label, parameter, map, dim);
        }

        @Override
        public double getDoubleValue() {
            if (dim == 0) {
                updateMap();
            }
            return parameter.getParameterValue(map[dim]);
        }

        private void updateMap() {
            Arrays.setAll(map, i -> i);
            Arrays.sort(map, (lhs, rhs) -> Double.compare(
                    parameter.getParameterValue(lhs), parameter.getParameterValue(rhs)));
        }
    }

    public enum Constraint {
        NONE,
        ORDERED
    }
}
