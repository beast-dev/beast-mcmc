package dr.evomodel.coalescent;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.AbstractModel;
import dr.util.Transform;

public class PopSizesLogger implements Loggable {
//    private final MultilocusNonparametricCoalescentLikelihood likelihood;
    private final Transform transform;
    private final String order;
    private final int popSizeDimension;
    private AbstractModel likelihood;

    public PopSizesLogger(MultilocusNonparametricCoalescentLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getPopSizeDimension();
        this.likelihood = likelihood;
    }

    // TODO this is a temporary constructor
    public PopSizesLogger(GMRFSkygridLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getDimension();
        this.likelihood = likelihood;
    }

    // TODO this is a temporary constructor
    public PopSizesLogger(GMRFMultilocusSkyrideLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getDimension();
        this.likelihood = likelihood;
    }

    private double getPopSize(int i) {
        if (likelihood instanceof MultilocusNonparametricCoalescentLikelihood) {
            return ((MultilocusNonparametricCoalescentLikelihood) likelihood).getPopulationSize(i);
        } else if (likelihood instanceof GMRFSkygridLikelihood) {
            double x = ((GMRFSkygridLikelihood) likelihood).getParameter().getParameterValue(i);
            return Math.exp(x);
        } else if (likelihood instanceof GMRFMultilocusSkyrideLikelihood) {
            double x = ((GMRFMultilocusSkyrideLikelihood) likelihood).getParameter().getParameterValue(i);
            return Math.exp(x);
        } else {
            throw new RuntimeException("Likelihood type not yet implemented.");
        }
    }
    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[popSizeDimension];
        for (int k = 0; k < popSizeDimension; ++k) {
            int finalK = k;
            columns[k] = new NumberColumn("PopSize" + "." + (finalK + 1)) {
                @Override
                public double getDoubleValue() {
                    if (transform != null) {
                        return transform.transform(getPopSize(finalK));
//                        return transform.transform(likelihood.getPopulationSize(finalK));
                    } else {
                        return getPopSize(finalK);
                    }
                }
            };
        }
        return columns;
    }
}
