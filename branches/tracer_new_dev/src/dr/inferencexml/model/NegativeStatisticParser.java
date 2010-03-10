package dr.inferencexml.model;

import dr.inference.model.NegativeStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class NegativeStatisticParser extends AbstractXMLObjectParser {

    public static String NEGATE_STATISTIC = "negativeStatistic";
    public static String NEGATE = "negate";
    public static String NEGATIVE = "negative";

    public String[] getParserNames() { return new String[] { getParserName(), NEGATIVE,  NEGATE}; }
    public String getParserName() { return NEGATE_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NegativeStatistic negativeStatistic = null;

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            negativeStatistic = new NegativeStatistic(NEGATE_STATISTIC, (Statistic)child);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return negativeStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise negation of the child statistic.";
    }

    public Class getReturnType() { return NegativeStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new ElementRule(Statistic.class, 1, 1 )
    };

}
