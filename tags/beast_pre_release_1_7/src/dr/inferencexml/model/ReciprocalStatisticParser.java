package dr.inferencexml.model;

import dr.inference.model.ReciprocalStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ReciprocalStatisticParser extends AbstractXMLObjectParser {

    public static String RECIPROCAL_STATISTIC = "reciprocalStatistic";
    public static String RECIPROCAL = "reciprocal";

    public String[] getParserNames() {
        return new String[]{getParserName(), RECIPROCAL};
    }

    public String getParserName() {
        return RECIPROCAL_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ReciprocalStatistic recipStatistic = null;

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            recipStatistic = new ReciprocalStatistic(xo.getId(), (Statistic) child);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return recipStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise reciprocal of the child statistic.";
    }

    public Class getReturnType() {
        return ReciprocalStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Statistic.class, 1, 1)
    };
}
