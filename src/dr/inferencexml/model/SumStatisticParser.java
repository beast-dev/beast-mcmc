package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.SumStatistic;
import dr.xml.*;

/**
 */
public class SumStatisticParser extends AbstractXMLObjectParser {

    public static String SUM_STATISTIC = "sumStatistic";
    public static String SUM = "sum";

    public String[] getParserNames() {
        return new String[]{getParserName(), SUM};
    }

    public String getParserName() {
        return SUM_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean elementwise = xo.getAttribute("elementwise", false);

        String name = SUM_STATISTIC;
        if (xo.hasAttribute(Statistic.NAME) || xo.hasAttribute(dr.xml.XMLParser.ID)) {
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        }

        final SumStatistic sumStatistic = new SumStatistic(name, elementwise);

        for (int i = 0; i < xo.getChildCount(); i++) {
            final Statistic statistic = (Statistic) xo.getChild(i);

            try {
                sumStatistic.addStatistic(statistic);
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Statistic added to " + getParserName() + " element is not of the same dimension");
            }
        }

        return sumStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise sum of the child statistics.";
    }

    public Class getReturnType() {
        return SumStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule("elementwise", true),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
    };
}
