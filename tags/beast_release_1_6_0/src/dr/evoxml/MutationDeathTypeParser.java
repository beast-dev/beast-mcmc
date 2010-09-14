package dr.evoxml;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.MutationDeathType;
import dr.evoxml.util.DataTypeUtils;
import dr.xml.*;

import java.util.logging.Logger;

/**
 *
 */
public class MutationDeathTypeParser extends AbstractXMLObjectParser {

    public static final String MODEL_NAME = "extendedDataType";
    public static final String STATE = "deathState";
    public static final String CODE = "code";
    public static final String STATES = "states";
    public static final String AMBIGUITY = "ambiguity";
    public static final String EXTANT = "extantState";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Logger.getLogger("dr.evolution").info("\nCreating an extended data type.");
        Logger.getLogger("dr.evolution").info("\tIf you publish results using this model, please reference Alekseyenko, Lee and Suchard (in submision).\n");

        DataType dataType = DataTypeUtils.getDataType(xo);

        char extantChar = '\0';
        XMLObject cxo = xo.getChild(EXTANT);
        if (cxo != null) {
            extantChar = cxo.getStringAttribute(CODE).charAt(0);
        }


        cxo = xo.getChild(STATE);
        char stateChar;
        if (cxo != null) {
            stateChar = cxo.getStringAttribute(CODE).charAt(0);
        } else {
            stateChar = dataType.getChar(dataType.getGapState());
        }

        Logger.getLogger("dr.evolution").info("\tNon-existent code: " + stateChar);

        if (dataType == null && extantChar == '\0')
            throw new XMLParseException("In " + xo.getName() + " you must either provide a data type or a code for extant state");

        MutationDeathType mdt;
        if (dataType != null) {
            Logger.getLogger("dr.evolution").info("\tBase type: " + dataType.toString());
            mdt = new MutationDeathType(dataType, stateChar);
        } else {
            Logger.getLogger("dr.evolution").info("\tExtant code: " + extantChar);
            mdt = new MutationDeathType(stateChar, extantChar);
        }

        char ambiguityChar = '\0';
        String states;
        Object dxo;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            dxo = xo.getChild(i);
            if (dxo instanceof XMLObject) {
                cxo = (XMLObject) dxo;
                if (cxo.getName().equals(AMBIGUITY)) {
                    ambiguityChar = cxo.getStringAttribute(CODE).charAt(0);
                    if (cxo.hasAttribute(STATES)) {
                        states = cxo.getStringAttribute(STATES);
                    } else {
                        states = "";
                    }
                    mdt.addAmbiguity(ambiguityChar, states);
                    Logger.getLogger("dr.evolution").info("\tAmbiguity code: " + ambiguityChar);
                }
            }
        }
        return mdt;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the MutationDeathType which extends a base datatype with an additional \"death\" state.";
    }

    public Class getReturnType() {
        return MutationDeathType.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(DataType.DATA_TYPE, "Base datatype name", DataType.getRegisteredDataTypeNames(), true),
            new ElementRule(DataType.class, true),
            new ElementRule(STATE, new XMLSyntaxRule[]{AttributeRule.newStringRule(CODE)}, true),
            new ElementRule(EXTANT, new XMLSyntaxRule[]{AttributeRule.newStringRule(CODE)}, true),
            new ContentRule("<ambiguity code=\"Z\" states=\"XY\"/>")
    };
}
