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
public class SequenceErrorModelComponent extends BaseComponentGenerator implements ComponentOptions {

    public SequenceErrorModelComponent(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {

        if (errorModelType == SequenceErrorType.NO_ERROR) {
            return false;
        }

        switch (point) {
            case AFTER_SITE_MODEL:
            case IN_TREE_LIKELIHOOD:
            case IN_FILE_LOG_PARAMETERS:
                return true;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        switch (point) {
            case AFTER_SITE_MODEL:
                writeErrorModel(writer);
                break;
            case IN_TREE_LIKELIHOOD:
                writer.writeTag(SequenceErrorModel.SEQUENCE_ERROR_MODEL,
                        new Attribute[]{new Attribute.Default<String>("idref", "errorModel")}, true);
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (errorModelType == SequenceErrorType.AGE_ALL || errorModelType == SequenceErrorType.AGE_TRANSITIONS) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "errorModel.ageRate"), true);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "errorModel.baseRate"), true);
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Sequence Error Model";
    }

    private void writeErrorModel(XMLWriter writer) {
        String errorType = (errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                errorModelType == SequenceErrorType.BASE_TRANSITIONS ?
                "transitions" : "all");


        writer.writeOpenTag(
                SequenceErrorModel.SEQUENCE_ERROR_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("id", "errorModel"),
                        new Attribute.Default<String>("type", errorType)
                }
        );

        if (errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                errorModelType == SequenceErrorType.AGE_ALL) {
            writeParameter(SequenceErrorModel.AGE_RELATED_RATE, "errorModel.ageRate", 1, writer);
        } else if (errorModelType == SequenceErrorType.BASE_TRANSITIONS ||
                errorModelType == SequenceErrorType.BASE_ALL) {
            writeParameter(SequenceErrorModel.BASE_ERROR_RATE, "errorModel.baseRate", 1, writer);
        } else {
            throw new IllegalArgumentException("Unknown error type");
        }

        writer.writeCloseTag(SequenceErrorModel.SEQUENCE_ERROR_MODEL);
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
        if (options.hasData()) {
            if (errorModelType == SequenceErrorType.AGE_ALL || errorModelType == SequenceErrorType.AGE_TRANSITIONS) {
                ops.add(modelOptions.getOperator("errorModel.ageRate"));
            } else if (errorModelType == SequenceErrorType.BASE_ALL || errorModelType == SequenceErrorType.BASE_TRANSITIONS) {
                ops.add(modelOptions.getOperator("errorModel.baseRate"));
            }
        }
    }

    public SequenceErrorType errorModelType = SequenceErrorType.NO_ERROR;
}
