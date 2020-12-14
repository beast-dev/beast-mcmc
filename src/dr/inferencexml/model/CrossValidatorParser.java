package dr.inferencexml.model;

import dr.inference.model.CrossValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class CrossValidatorParser extends AbstractXMLObjectParser {

    public final static String LOG_SUM = "logSum";
    private final static String TYPE = "type";
    public final static String CROSS_VALIDATION = "crossValidation";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CrossValidationProvider provider = (CrossValidationProvider) xo.getChild(CrossValidationProvider.class);
        boolean logSum = xo.getAttribute(LOG_SUM, false);

        CrossValidationProvider.ValidationType validationType;

        String validation = xo.getAttribute(TYPE, CrossValidationProvider.ValidationType.SQUARED_ERROR.getName());

        if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.SQUARED_ERROR.getName())) {
            validationType = CrossValidationProvider.ValidationType.SQUARED_ERROR;
        } else if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.BIAS.getName())) {
            validationType = CrossValidationProvider.ValidationType.BIAS;
        } else if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.VALUE.getName())) {
            validationType = CrossValidationProvider.ValidationType.VALUE;
        } else {
            throw new XMLParseException("The attribute '" + TYPE + "' can only take values '" +
                    CrossValidationProvider.ValidationType.SQUARED_ERROR.getName() + "', " +
                    CrossValidationProvider.ValidationType.BIAS.getName() + "', or" +
                    CrossValidationProvider.ValidationType.VALUE.getName() + "'.");
        }

        if (logSum) return new CrossValidationProvider.CrossValidatorSum(provider, validationType);
        return new CrossValidationProvider.CrossValidator(provider, validationType);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(LOG_SUM, true),
                AttributeRule.newStringRule(TYPE, true),
                new ElementRule(CrossValidationProvider.class)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CrossValidationProvider.CrossValidator.class;
    }

    @Override
    public String getParserName() {
        return CROSS_VALIDATION;
    }
}
