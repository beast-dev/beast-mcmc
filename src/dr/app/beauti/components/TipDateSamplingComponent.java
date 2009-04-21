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
public class TipDateSamplingComponent extends BaseComponentGenerator implements ComponentOptions {

    public TipDateSamplingComponent(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {

        if (tipDateSamplingType == TipDateSamplingType.NO_SAMPLING) {
            return false;
        }

        switch (point) {
            case IN_TREE_MODEL:
            case IN_MCMC_PRIOR:
            case IN_FILE_LOG_PARAMETERS:
                return true;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        switch (point) {
            case IN_TREE_MODEL:
                break;
            case IN_MCMC_PRIOR:
                break;
            case IN_FILE_LOG_PARAMETERS:
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "treeModel.tipDates"), true);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Tip date sampling";
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter("treeModel.tipDates", "date of specified tips", ModelOptions.TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        modelOptions.createScaleOperator("treeModel.tipDates", 3.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_ALL) {
            params.add(modelOptions.getParameter("errorModel.ageRate"));
        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_SET) {
            params.add(modelOptions.getParameter("errorModel.baseRate"));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_ALL) {
            ops.add(modelOptions.getOperator("errorModel.ageRate"));
        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_SET) {
            ops.add(modelOptions.getOperator("errorModel.baseRate"));
        }
    }

    public TipDateSamplingType tipDateSamplingType = TipDateSamplingType.NO_SAMPLING;
}