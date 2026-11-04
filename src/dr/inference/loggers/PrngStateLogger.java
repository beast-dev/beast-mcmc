package dr.inference.loggers;

import dr.math.MathUtils;

public class PrngStateLogger implements Loggable {

    private static final String LABEL = "prng";

    @Override
    public LogColumn[] getColumns() {
        int dim = MathUtils.getRandomState().length;
        LogColumn[] columns = new LogColumn[dim];

        for (int i = 0; i < dim; ++i) {
            final int index = i;
            columns[i] = new LogColumn.Abstract(LABEL + (index + 1)) {

                @Override
                protected String getFormattedValue() {
                    return Integer.toString(MathUtils.getRandomState()[index]);
                }
            };
        }

        return columns;
    }
}
