package dr.app.beauti.components.sequenceerror;

import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.SequenceErrorType;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentOptions implements ComponentOptions {
    static public final String ERROR_MODEL = "errorModel";
    static public final String AGE_RATE = "ageRelatedErrorRate";
    static public final String BASE_RATE = "baseErrorRate";

    static public final String AGE_RATE_PARAMETER = ERROR_MODEL + "." + AGE_RATE;
    static public final String BASE_RATE_PARAMETER = ERROR_MODEL + "." + BASE_RATE;

    static public final String HYPERMUTION_RATE_PARAMETER = "hypermutation.rate";
    static public final String HYPERMUTANT_INDICATOR_PARAMETER = "hypermutant.indicator";
    static public final String HYPERMUTANT_COUNT_STATISTIC = "hypermutation.count";

    public SequenceErrorModelComponentOptions() {
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createNonNegativeParameterInfinitePrior(AGE_RATE_PARAMETER,"age dependent sequence error rate",
                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0E-8);
        modelOptions.createZeroOneParameterUniformPrior(BASE_RATE_PARAMETER,"base sequence error rate", 1.0E-8);

        modelOptions.createZeroOneParameterUniformPrior(HYPERMUTION_RATE_PARAMETER,"APOBEC editing rate per context", 1.0E-8);
        modelOptions.createParameter(HYPERMUTANT_INDICATOR_PARAMETER, "indicator parameter reflecting which sequences are hypermutated", 0.0);

        modelOptions.createDiscreteStatistic(HYPERMUTANT_COUNT_STATISTIC, "count of the number of hypermutated sequences");

        modelOptions.createScaleOperator(AGE_RATE_PARAMETER, modelOptions.demoTuning, 3.0);
        modelOptions.createOperator(BASE_RATE_PARAMETER, OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);

        modelOptions.createOperator(HYPERMUTION_RATE_PARAMETER, OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);
        modelOptions.createOperator(HYPERMUTANT_INDICATOR_PARAMETER, OperatorType.BITFLIP, -1.0, 10);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        if (isHypermutation()) {
            params.add(modelOptions.getParameter(HYPERMUTION_RATE_PARAMETER));
            params.add(modelOptions.getParameter(HYPERMUTANT_INDICATOR_PARAMETER));
            params.add(modelOptions.getParameter(HYPERMUTANT_COUNT_STATISTIC));
        }
        if (hasAgeDependentRate()) {
            params.add(modelOptions.getParameter(AGE_RATE_PARAMETER));
        }
        if (hasBaseRate()) {
            params.add(modelOptions.getParameter(BASE_RATE_PARAMETER));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        if (isHypermutation()) {
            ops.add(modelOptions.getOperator(HYPERMUTION_RATE_PARAMETER));
            ops.add(modelOptions.getOperator(HYPERMUTANT_INDICATOR_PARAMETER));
        }
        if (hasAgeDependentRate()) {
            ops.add(modelOptions.getOperator(AGE_RATE_PARAMETER));
        }
        if (hasBaseRate()) {
            ops.add(modelOptions.getOperator(BASE_RATE_PARAMETER));
        }
    }

    public boolean hasAgeDependentRate() {
        return (errorModelType == SequenceErrorType.AGE_ALL) || (errorModelType == SequenceErrorType.AGE_TRANSITIONS);
    }

    public boolean hasBaseRate() {
        return (errorModelType == SequenceErrorType.BASE_ALL) || (errorModelType == SequenceErrorType.BASE_TRANSITIONS);
    }

    public boolean isHypermutation() {
        return (errorModelType == SequenceErrorType.HYPERMUTATION_ALL) ||
                (errorModelType == SequenceErrorType.HYPERMUTATION_BOTH) ||
        (errorModelType == SequenceErrorType.HYPERMUTATION_HA3G) ||
        (errorModelType == SequenceErrorType.HYPERMUTATION_HA3F);
    }

    public SequenceErrorType errorModelType = SequenceErrorType.NO_ERROR;
}