package dr.evomodel.speciation;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Statistic;

public class BirthDeathCompoundParameterLogger extends Statistic.Abstract {

    private final NewBirthDeathSerialSamplingModel bdss;
    private final int dim;

    public BirthDeathCompoundParameterLogger(NewBirthDeathSerialSamplingModel bdss) {
        this.bdss = bdss;
        this.dim = bdss.getDeathRateParameter().getDimension();
    }

    private double getCompoundParameter(int i) {
        double birth = bdss.getBirthRateParameter().getParameterValue(i);
        double death = bdss.getDeathRateParameter().getParameterValue(i);
        double sampling = bdss.getSamplingRateParameter().getParameterValue(i);
        double treatment = bdss.getTreatmentProbabilityParameter().getParameterValue(i);

        return birth / (sampling * treatment + death);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double getStatisticValue(int dim) {
        return getCompoundParameter(dim);
    }
}
