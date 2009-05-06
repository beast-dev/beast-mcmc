package dr.app.beauti.components;

import dr.app.beauti.options.*;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentOptions implements ComponentOptions {
    static public final String ERROR_MODEL = "errorModel";
    static private final String AGE_RATE = "ageRate";
    static private final String BASE_RATE = "baseRate";

    public SequenceErrorModelComponentOptions() {
    }

    public String errTypeName() {
        switch( errorModelType ) {
            case AGE_ALL:
            case AGE_TRANSITIONS:
                return AGE_RATE;
            case  BASE_ALL:
            case BASE_TRANSITIONS:
                return BASE_RATE;
        }
        return null;
    }

    public String qualifiedErrTypeName() {
        return ERROR_MODEL + "." + errTypeName();
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter(ERROR_MODEL + "." + AGE_RATE,
                "age dependent sequence error rate", ModelOptions.SUBSTITUTION_RATE_SCALE, 1.0E-8, 0.0, Double.POSITIVE_INFINITY);
        modelOptions.createParameter(ERROR_MODEL + "." + BASE_RATE,
                "base sequence error rate", ModelOptions.UNITY_SCALE, 1.0E-8, 0.0, 1.0);

        modelOptions.createScaleOperator(ERROR_MODEL + "." + AGE_RATE, 3.0);
        modelOptions.createOperator(ERROR_MODEL + "." + BASE_RATE, OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        params.add(modelOptions.getParameter(qualifiedErrTypeName()));
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        ops.add(modelOptions.getOperator(qualifiedErrTypeName()));
    }

    public static SequenceErrorType errorModelType = SequenceErrorType.NO_ERROR;
}