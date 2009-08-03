package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.evomodel.treelikelihood.SequenceErrorModel;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentGenerator extends BaseComponentGenerator {

    public SequenceErrorModelComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);

        if (comp.errorModelType == SequenceErrorType.NO_ERROR) {
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
        SequenceErrorModelComponentOptions component = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);

        switch (point) {
            case AFTER_SITE_MODEL:
                writeErrorModel(writer, component);
                break;
            case IN_TREE_LIKELIHOOD:
            	writer.writeIDref(SequenceErrorModel.SEQUENCE_ERROR_MODEL, "errorModel");
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (component.hasAgeDependentRate()) {
                	writer.writeIDref(ParameterParser.PARAMETER, SequenceErrorModelComponentOptions.AGE_RATE_PARAMETER);
                }
                if (component.hasBaseRate()) {
                	writer.writeIDref(ParameterParser.PARAMETER, SequenceErrorModelComponentOptions.BASE_RATE_PARAMETER);
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Sequence Error Model";
    }

    private void writeErrorModel(XMLWriter writer, SequenceErrorModelComponentOptions component) {
        final String errorType = (component.errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                component.errorModelType == SequenceErrorType.BASE_TRANSITIONS ?
                "transitions" : "all");


        writer.writeOpenTag(
                SequenceErrorModel.SEQUENCE_ERROR_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, SequenceErrorModelComponentOptions.ERROR_MODEL),
                        new Attribute.Default<String>("type", errorType)
                }
        );

        if (component.hasAgeDependentRate()) {
            writeParameter(SequenceErrorModelComponentOptions.AGE_RATE, SequenceErrorModelComponentOptions.AGE_RATE_PARAMETER, 1, writer);
        }
        if (component.hasBaseRate()) {
            writeParameter(SequenceErrorModelComponentOptions.BASE_RATE, SequenceErrorModelComponentOptions.BASE_RATE_PARAMETER, 1, writer);
        }

        writer.writeCloseTag(SequenceErrorModel.SEQUENCE_ERROR_MODEL);
    }
}
