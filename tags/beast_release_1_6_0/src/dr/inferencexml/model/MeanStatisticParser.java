package dr.inferencexml.model;

import dr.inference.model.MeanStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class MeanStatisticParser extends AbstractXMLObjectParser {

    public static String MEAN_STATISTIC = "meanStatistic";

    public String getParserName() { return MEAN_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MeanStatistic meanStatistic = new MeanStatistic(MEAN_STATISTIC);

        for (int i =0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Statistic) {
                meanStatistic.addStatistic((Statistic)child);
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        return meanStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the mean of the child statistics.";
    }

    public Class getReturnType() { return MeanStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
    };

}
