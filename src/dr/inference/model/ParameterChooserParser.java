package dr.inference.model;

import dr.xml.*;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ParameterChooserParser extends dr.xml.AbstractXMLObjectParser {
    public static String VARIABLE_SELECTOR = "variableSelector";
    public static String INDEX = "index";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String name = xo.getId();
        final ValuesPool pool = (ValuesPool) xo.getChild(ValuesPool.class);
        final int[] which = xo.getIntegerArrayAttribute(INDEX);
        for( int w : which ) {
            if( ! (0 <= w && w < pool.length()) ) {
                throw new XMLParseException("index " + w + " out of range");
            }
        }
        return new ParameterChooser(name, pool, which);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(ValuesPool.class,1,1),
                AttributeRule.newIntegerArrayRule(INDEX, false),
        };
    }

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return ParameterChooser.class;
    }

    public String getParserName() {
        return VARIABLE_SELECTOR;
    }
}
