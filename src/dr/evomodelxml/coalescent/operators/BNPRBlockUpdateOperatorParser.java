package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.operators.BNPRBlockUpdateOperator;
import dr.xml.AbstractXMLObjectParser;

/**
 * Created by mkarcher on 11/3/16.
 */
public class BNPRBlockUpdateOperatorParser extends AbstractXMLObjectParser {
    public static final String BNPR_BLOCK_OPERATOR = "bnprBlockUpdateOperator";

    public String getParserName() {
        return BNPR_BLOCK_OPERATOR;
    }

    public Class getReturnType() {
        return BNPRBlockUpdateOperator.class;
    }
}
