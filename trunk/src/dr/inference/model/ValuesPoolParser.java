package dr.inference.model;

import dr.xml.*;
//import sun.awt.X11.XLayerProtocol;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ValuesPoolParser extends dr.xml.AbstractXMLObjectParser {
    public static String VALUES_POOL = "valuesPool";
    public static String VALUES = "values";
     public static String SELECTOR = "selector";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Variable pool = (Variable)xo.getElementFirstChild(VALUES);
        final Variable selector = (Variable)xo.getElementFirstChild(SELECTOR);

        if( pool.getSize() != selector.getSize() ) {
            throw new XMLParseException("variable Pool and selector must have equal length.");
        }
        return new ValuesPool(pool, selector);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
          return new XMLSyntaxRule[] {
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
