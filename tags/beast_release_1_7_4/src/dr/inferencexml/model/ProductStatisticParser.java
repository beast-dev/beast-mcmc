package dr.inferencexml.model;

import dr.inference.model.ProductStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ProductStatisticParser extends AbstractXMLObjectParser {

    public static String PRODUCT_STATISTIC = "productStatistic";
    public static String PRODUCT = "product";
    public static String ELEMENT_WISE = "elementwise";

    public String[] getParserNames() {
        return new String[]{getParserName(), PRODUCT};
    }

    public String getParserName() {
        return PRODUCT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        boolean elementwise = xo.getAttribute(ELEMENT_WISE, false);

        ProductStatistic productStatistic = new ProductStatistic(PRODUCT_STATISTIC, elementwise);

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Statistic) {
                try {
                    productStatistic.addStatistic((Statistic) child);
                } catch (IllegalArgumentException iae) {
                    throw new XMLParseException("Statistic added to " + getParserName() + " element is not of the same dimension");
                }
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        return productStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the product of the child statistics.";
    }

    public Class getReturnType() {
        return ProductStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ELEMENT_WISE, true),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
    };
}
