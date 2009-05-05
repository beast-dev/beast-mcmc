package dr.app.beauti.components;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.evomodel.treelikelihood.SequenceErrorModel;
import dr.util.Attribute;
import dr.inference.model.ParameterParser;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentOptions implements ComponentOptions {

    public SequenceErrorModelComponentOptions() {
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter("errorModel.ageRate", "age dependent sequence error rate", ModelOptions.SUBSTITUTION_RATE_SCALE, 1.0E-8, 0.0, Double.POSITIVE_INFINITY);
        modelOptions.createParameter("errorModel.baseRate", "base sequence error rate", ModelOptions.UNITY_SCALE, 1.0E-8, 0.0, 1.0);

        modelOptions.createScaleOperator("errorModel.ageRate", 3.0);
        modelOptions.createOperator("errorModel.baseRate", OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        if (errorModelType == SequenceErrorType.AGE_ALL || errorModelType == SequenceErrorType.AGE_TRANSITIONS) {
            params.add(modelOptions.getParameter("errorModel.ageRate"));
        } else if (errorModelType == SequenceErrorType.BASE_ALL || errorModelType == SequenceErrorType.BASE_TRANSITIONS) {
            params.add(modelOptions.getParameter("errorModel.baseRate"));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        if (errorModelType == SequenceErrorType.AGE_ALL || errorModelType == SequenceErrorType.AGE_TRANSITIONS) {
            ops.add(modelOptions.getOperator("errorModel.ageRate"));
        } else if (errorModelType == SequenceErrorType.BASE_ALL || errorModelType == SequenceErrorType.BASE_TRANSITIONS) {
            ops.add(modelOptions.getOperator("errorModel.baseRate"));
        }
    }

    public SequenceErrorType errorModelType = SequenceErrorType.NO_ERROR;
}