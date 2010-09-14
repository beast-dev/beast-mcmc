package dr.inferencexml.model;

import dr.inference.model.ExponentialStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ExponentialStatisticParser extends AbstractXMLObjectParser {

    public static String EXPONENTIAL_STATISTIC = "exponentialStatistic";
    public static String EXP = "exp";

    public String[] getParserNames() {
        return new String[]{getParserName(), EXP};
    }

    public String getParserName() {
        return EXPONENTIAL_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ExponentialStatistic expStatistic = null;

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            expStatistic = new ExponentialStatistic(EXPONENTIAL_STATISTIC, (Statistic) child);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return expStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise exponentiation of the child statistic.";
    }

    public Class getReturnType() {
        return ExponentialStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Statistic.class, 1, 1)
    };
}
