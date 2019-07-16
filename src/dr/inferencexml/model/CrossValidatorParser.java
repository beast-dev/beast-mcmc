package dr.inferencexml.model;

import dr.inference.model.CrossValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class CrossValidatorParser extends AbstractXMLObjectParser {

    final static String LOG_SUM = "logSum";
    final static String CROSS_VALIDATION = "crossValidation";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CrossValidationProvider provider = (CrossValidationProvider) xo.getChild(CrossValidationProvider.class);
        boolean logSum = xo.getAttribute(LOG_SUM, false);

        if (logSum) return new CrossValidationProvider.CrossValidatorSum(provider);
        return new CrossValidationProvider.CrossValidator(provider);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(LOG_SUM, true),
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
