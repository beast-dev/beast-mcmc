package dr.evomodel.speciation;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Statistic;
import dr.math.UnivariateMinimum;

// TODO use "implements Loggable" instead of "extends Statistic.Abstract"
// TODO hard-code names
public class BirthDeathCompoundParameterLogger extends Statistic.Abstract {

    public enum BDPCompoundParameterType {
        EFFECTIVE_REPRODUCTIVE_NUMBER("effectiveReproductiveNumber") {
            public double getCompoundParameterForType(double birthRate, double deathRate, double samplingRate, double treatmentProbability, double samplingProbability) {
                return birthRate / (samplingRate * treatmentProbability + deathRate);
            }
        };

        BDPCompoundParameterType(String name) {
            this.name = name;
//            this.label = label;
        }

        public String getName() {
            return name;
        }

//        public String getLabel() {
//            return label;
//        }

        // TODO may need to be
        public abstract double getCompoundParameterForType(double birthRate, double deathRate, double samplingRate, double treatmentProbability, double samplingProbability);

        private String name;
//        private String label;
    }

    public BirthDeathCompoundParameterLogger(NewBirthDeathSerialSamplingModel bdss, BDPCompoundParameterType type) {
        this.bdss = bdss;
        this.type = type;
        this.dim = bdss.getDeathRateParameter().getDimension();
    }

    private double getCompoundParameter(int i) {
        double birth = bdss.getBirthRateParameter().getParameterValue(i);
        double death = bdss.getDeathRateParameter().getParameterValue(i);
        double sampling = bdss.getSamplingRateParameter().getParameterValue(i);
        double treatment = bdss.getTreatmentProbabilityParameter().getParameterValue(i);
        double prob = bdss.getSamplingProbabilityParameter().getParameterValue(i);

        return type.getCompoundParameterForType(birth, death, sampling, treatment, prob);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double getStatisticValue(int dim) {
        return getCompoundParameter(dim);
    }

    private final NewBirthDeathSerialSamplingModel bdss;
    private final int dim;
    private final BDPCompoundParameterType type;

}
