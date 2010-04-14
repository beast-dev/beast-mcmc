package dr.app.beauti.components;

import dr.app.beauti.enumTypes.SequenceErrorType;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.XMLWriter;
import dr.evomodelxml.treelikelihood.SequenceErrorModelParser;
import dr.inference.model.*;
import dr.util.Attribute;
import dr.xml.XMLParser;
import dr.evomodel.treelikelihood.APOBECErrorModel;
import dr.inferencexml.model.SumStatisticParser;

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
                writer.writeIDref(SequenceErrorModelParser.SEQUENCE_ERROR_MODEL, "errorModel");
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (component.isHypermutation()) {
                    writer.writeIDref(ParameterParser.PARAMETER, SequenceErrorModelComponentOptions.HYPERMUTION_RATE_PARAMETER);
                    writer.writeIDref(StatisticParser.STATISTIC, SequenceErrorModelComponentOptions.HYPERMUTANT_COUNT_STATISTIC);
                    writer.writeOpenTag(StatisticParser.STATISTIC,
                            new Attribute.Default<String>("name", "isHypermutated"));
                    writer.writeIDref(APOBECErrorModel.APOBEC_ERROR_MODEL,
                            SequenceErrorModelComponentOptions.ERROR_MODEL);
                    writer.writeCloseTag(StatisticParser.STATISTIC);
                }

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
        if (component.isHypermutation()) {

            final String errorType;
            switch (component.errorModelType) {
                case HYPERMUTATION_ALL:
                    errorType = APOBECErrorModel.APOBECType.ALL.toString();
                    break;
                case HYPERMUTATION_BOTH:
                    errorType = APOBECErrorModel.APOBECType.BOTH.toString();
                    break;
                case HYPERMUTATION_HA3F:
                    errorType = APOBECErrorModel.APOBECType.HA3G.toString();
                    break;
                case HYPERMUTATION_HA3G:
                    errorType = APOBECErrorModel.APOBECType.HA3F.toString();
                    break;
                default:
                    throw new RuntimeException("Unknown ErrorModelType: " + component.errorModelType.toString());
            }
            writer.writeOpenTag(
                    APOBECErrorModel.APOBEC_ERROR_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, SequenceErrorModelComponentOptions.ERROR_MODEL),
                            new Attribute.Default<String>("type", errorType)
                    }
            );

            writer.writeIDref("taxa", "taxa");


            writeParameter(APOBECErrorModel.HYPERMUTATION_RATE, SequenceErrorModelComponentOptions.HYPERMUTION_RATE_PARAMETER, 1, writer);
            writeParameter(APOBECErrorModel.HYPERMUTATION_INDICATORS, SequenceErrorModelComponentOptions.HYPERMUTANT_INDICATOR_PARAMETER, 1, writer);

            writer.writeCloseTag(APOBECErrorModel.APOBEC_ERROR_MODEL);

            writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, SequenceErrorModelComponentOptions.HYPERMUTANT_COUNT_STATISTIC),
                    new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
            writer.writeIDref(ParameterParser.PARAMETER, SequenceErrorModelComponentOptions.HYPERMUTANT_INDICATOR_PARAMETER);
            writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);
        } else {
            final String errorType = (component.errorModelType == SequenceErrorType.AGE_TRANSITIONS ||
                    component.errorModelType == SequenceErrorType.BASE_TRANSITIONS ?
                    "transitions" : "all");


            writer.writeOpenTag(
                    SequenceErrorModelParser.SEQUENCE_ERROR_MODEL,
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

            writer.writeCloseTag(SequenceErrorModelParser.SEQUENCE_ERROR_MODEL);
        }
    }
}
