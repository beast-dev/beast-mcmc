package dr.inferencexml.model;

import dr.inference.model.NotStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class NotStatisticParser extends AbstractXMLObjectParser {

    public static String NOT_STATISTIC = "notStatistic";
    public static String NOT = "not";

    public String[] getParserNames() {
        return new String[]{getParserName(), NOT};
    }

    public String getParserName() {
        return NOT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NotStatistic notStatistic;

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            notStatistic = new NotStatistic(NOT_STATISTIC, (Statistic) child);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return notStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise inverse of the child statistic.";
    }

    public Class getReturnType() {
        return NotStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Statistic.class, 1, 1)
    };
}
