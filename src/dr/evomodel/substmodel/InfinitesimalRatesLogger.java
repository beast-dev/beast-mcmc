package dr.evomodel.substmodel;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;

public class InfinitesimalRatesLogger implements Loggable {

    public InfinitesimalRatesLogger(SubstitutionModel substitutionModel) {
        this.substitutionModel = substitutionModel;
    }

    @Override
    public LogColumn[] getColumns() {
        int stateCount = substitutionModel.getDataType().getStateCount();

        if (generator == null) {
            generator = new double[stateCount * stateCount];
        }

        LogColumn[] columns = new LogColumn[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                final int k = i * stateCount + j;
                columns[k] = new NumberColumn(substitutionModel.getId() + "." + (i + 1) + "." + (j + 1)) {
                    @Override
                    public double getDoubleValue() {
                        if (k == 0) { // Refresh at first-element read
                            substitutionModel.getInfinitesimalMatrix(generator);
                        }
                        return generator[k];
                    }
                };
            }
        }

        return columns;
    }

    private final SubstitutionModel substitutionModel;
    private double[] generator;
}
