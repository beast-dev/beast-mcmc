package dr.inferencexml.model;

import dr.inference.model.IndicatorProductParameter;
import dr.inference.model.Parameter;

import dr.xml.*;


public class IndicatorProductParameterParser extends AbstractXMLObjectParser {

    public static final String INDICATOR_PRODUCT_PARAMETER = "indicatorProductParameter";
    public static final String INDICATOR = "indicator";
    public static final String CONTINUOUS_PARAMETER = "continuousParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter continuous = (Parameter) xo.getChild(CONTINUOUS_PARAMETER).getChild(Parameter.class);
        Parameter indicator = (Parameter) xo.getChild(INDICATOR).getChild(Parameter.class);

        IndicatorProductParameter prodParam = new IndicatorProductParameter(indicator, continuous);

        return prodParam;


    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(INDICATOR, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, false),
        new ElementRule(CONTINUOUS_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, false)
    };

    public String getParserDescription() {
        return "A element-wise product of an indicator and continuous paramter.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return INDICATOR_PRODUCT_PARAMETER;
    }
}
