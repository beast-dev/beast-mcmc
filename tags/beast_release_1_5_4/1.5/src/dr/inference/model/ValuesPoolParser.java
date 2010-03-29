package dr.inference.model;

import dr.xml.*;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ValuesPoolParser extends dr.xml.AbstractXMLObjectParser {
    public static String VALUES_POOL = "valuesPool";
    public static String VALUES = "values";
    public static String SELECTOR = "selector";
    private static final String DEFAULT_VALUE = "default";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Variable<Double> pool = (Variable<Double>)xo.getElementFirstChild(VALUES);
        final Variable<Double> selector = (Variable<Double>)xo.getElementFirstChild(SELECTOR);
        final double d = xo.getDoubleAttribute(DEFAULT_VALUE);

        if( pool.getSize() != selector.getSize() ) {
            throw new XMLParseException("variable Pool and selector must have equal length.");
        }
        return new ValuesPool(pool, selector, d);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(DEFAULT_VALUE),
                new ElementRule(VALUES, new XMLSyntaxRule[]{
                        new ElementRule(Variable.class,1,1) }),
                new ElementRule(SELECTOR, new XMLSyntaxRule[]{
                        new ElementRule(Variable.class,1,1) })
          };
      }


    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return ValuesPool.class;
    }

    public String getParserName() {
        return VALUES_POOL;
    }
}
