package dr.app.beauti.components;

import dr.app.beauti.XMLWriter;
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
        SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions)options.getComponentOptions(SequenceErrorModelComponentOptions.class);

        switch (point) {
            case AFTER_SITE_MODEL:
                writeErrorModel(writer, comp);
                break;
            case IN_TREE_LIKELIHOOD:
                writer.writeTag(SequenceErrorModel.SEQUENCE_ERROR_MODEL,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "errorModel")}, true);
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (comp.errorModelType == SequenceErrorType.AGE_ALL || comp.errorModelType == SequenceErrorType.AGE_TRANSITIONS) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "errorModel.ageRate"), true);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "errorModel.baseRate"), true);
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
        String errorType = (component.errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                component.errorModelType == SequenceErrorType.BASE_TRANSITIONS ?
                "transitions" : "all");


        writer.writeOpenTag(
                SequenceErrorModel.SEQUENCE_ERROR_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "errorModel"),
                        new Attribute.Default<String>("type", errorType)
                }
        );

        if (component.errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                component.errorModelType == SequenceErrorType.AGE_ALL) {
            writeParameter(SequenceErrorModel.AGE_RELATED_RATE, "errorModel.ageRate", 1, writer);
        } else if (component.errorModelType == SequenceErrorType.BASE_TRANSITIONS ||
                component.errorModelType == SequenceErrorType.BASE_ALL) {
            writeParameter(SequenceErrorModel.BASE_ERROR_RATE, "errorModel.baseRate", 1, writer);
        } else {
            throw new IllegalArgumentException("Unknown error type");
        }

        writer.writeCloseTag(SequenceErrorModel.SEQUENCE_ERROR_MODEL);
    }

}
