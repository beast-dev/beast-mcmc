package dr.evoxml;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Package: KStateDataTypeParser
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Oct 13, 2009
 *         Time: 3:42:14 PM
 */
public class KStateDataTypeParser extends AbstractXMLObjectParser {
    public static final String K_STATE_DATATYPE = "kStateType";
    public static final String STATE_COUNT = "stateCount";
    public static final String START_WITH = "startWith";

    //public static XMLObjectParser PARSER=new KStateDataTypeParser();

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int k = xo.getIntegerAttribute(STATE_COUNT);
        int sw = 0;
        if (xo.hasAttribute(START_WITH))
            sw = xo.getIntegerAttribute(START_WITH);
        Collection<String> states = new ArrayList<String>();

        System.err.println(states.toArray().toString());

        for (int i = sw; i < k + sw; ++i) {
            states.add(Integer.toString(i));
        }
        System.err.println("Constructing " + k + "-state datatype");
        return new GeneralDataType(states);
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(STATE_COUNT, false),
            AttributeRule.newIntegerRule(START_WITH, true)};

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "Parser for k-state model.";
    }

    public Class getReturnType() {
        return DataType.class;
    }

    public String getParserName() {
        return K_STATE_DATATYPE;
    }

}
