package dr.inferencexml.model;

import dr.inference.model.RatioStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class RatioStatisticParser extends AbstractXMLObjectParser {

    public static String RATIO_STATISTIC = "ratioStatistic";
    public static String RATIO = "ratio";

    public String[] getParserNames() { return new String[] { getParserName(), RATIO}; }
    public String getParserName() { return RATIO_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Statistic numerator = (Statistic)xo.getChild(0);
        Statistic denominator = (Statistic)xo.getChild(1);

        RatioStatistic ratioStatistic;
        try {
            ratioStatistic = new RatioStatistic(RATIO_STATISTIC, numerator, denominator);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("Error parsing " + getParserName() + " the numerator and denominator statistics " +
                    "should be of the same dimension or of dimension 1");
        }
        return ratioStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the ratio of the 2 child statistics.";
    }

    public Class getReturnType() { return RatioStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Statistic.class, "The two operand statistics", 2, 2)
    };

}
