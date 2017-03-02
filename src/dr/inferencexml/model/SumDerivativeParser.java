package dr.inferencexml.model;

import dr.inference.model.PotentialDerivativeInterface;
import dr.inference.model.SumDerivative;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Max Tolkoff
 */
public class SumDerivativeParser extends AbstractXMLObjectParser{
    public final static String SUM_DERIVATIVE = "sumDerivative";

    @Override
    public String getParserName() {
        return SUM_DERIVATIVE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ArrayList<PotentialDerivativeInterface> derivativeList = new ArrayList<PotentialDerivativeInterface>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            derivativeList.add((PotentialDerivativeInterface) xo.getChild(i));
        }


        return new SumDerivative(derivativeList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PotentialDerivativeInterface.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SumDerivative.class;
    }
}
