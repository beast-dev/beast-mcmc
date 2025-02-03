package dr.evomodel.coalescent;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.util.Transform;

public class PopSizesLogger implements Loggable {
    private final MultilocusNonparametricCoalescentLikelihood likelihood;
    private final Transform transform;
    private final String order;

    public PopSizesLogger(MultilocusNonparametricCoalescentLikelihood likelihood, Transform transform, String order) {
        this.likelihood = likelihood;
        this.transform = transform;
        this.order = order; // not yet implemented
    }

    @Override
    public LogColumn[] getColumns() {
        int nOutputs = likelihood.getNGridPoints() + 1;

        LogColumn[] columns = new LogColumn[nOutputs];
        for (int k = 0; k < nOutputs; ++k) {
            int finalK = k;
            columns[k] = new NumberColumn("PopSize" + "." + (finalK + 1)) {
                @Override
                public double getDoubleValue() {
                    if (transform != null) {
                        return transform.transform(likelihood.getPopulationSize(finalK));
                    } else {
                        return likelihood.getPopulationSize(finalK);
                    }
                }
            };
        }
        return columns;
    }
}
