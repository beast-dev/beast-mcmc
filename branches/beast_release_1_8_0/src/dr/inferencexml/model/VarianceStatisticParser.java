package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.VarianceStatistic;
import dr.xml.*;

/**
 */
public class VarianceStatisticParser extends AbstractXMLObjectParser {

    public static String VARIANCE_STATISTIC = "varianceStatistic";

    public String getParserName() { return VARIANCE_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        VarianceStatistic varStatistic = new VarianceStatistic(VARIANCE_STATISTIC);

        for (int i =0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Statistic) {
                varStatistic.addStatistic((Statistic)child);
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        return varStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the variance of the child statistics.";
    }

    public Class getReturnType() { return VarianceStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
    };

}
