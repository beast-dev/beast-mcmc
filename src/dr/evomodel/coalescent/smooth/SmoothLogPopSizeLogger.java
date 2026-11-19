package dr.evomodel.coalescent.smooth;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;


/**
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SmoothLogPopSizeLogger implements Loggable {

    private final NewSmoothSkygridLikelihood likelihood;

    public SmoothLogPopSizeLogger(NewSmoothSkygridLikelihood likelihood) {
        this.likelihood = likelihood;
    }

    @Override
    public LogColumn[] getColumns() {
        int size = likelihood.getGridPointParameter().getDimension();

        LogColumn[] array = new LogColumn[size];

        for (int i = 0; i < size; ++i) {

            final int index = i;
            array[i] = new NumberColumn("log.smooth.pop.size." + (index + 1)) {
                @Override
                public double getDoubleValue() {
                    return -Math.log(likelihood.getSmoothPopSizeInverse(likelihood.getGridPointParameter().getParameterValue(index)));
                }
            };
        }

        return array;
    }
}
