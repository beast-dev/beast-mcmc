package dr.inference.operators;

import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * @author Alexei Drummond
 */
public enum CoercionMode {
    DEFAULT, COERCION_ON, COERCION_OFF;

    public static CoercionMode parseMode(XMLObject xo) throws XMLParseException {
        CoercionMode mode = CoercionMode.DEFAULT;
        if (xo.hasAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
            if (xo.getBooleanAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                mode = CoercionMode.COERCION_ON;
            } else {
                mode = CoercionMode.COERCION_OFF;
            }
        }
        return mode;
    }
}
