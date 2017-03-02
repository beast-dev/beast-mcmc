package dr.inferencexml.model;

import dr.inference.model.AppendedPotentialDerivative;
import dr.inference.model.PotentialDerivativeInterface;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivativeParser extends AbstractXMLObjectParser {

    public final static String SUM_DERIVATIVE = "appendedPotentialDerivative";

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


        return new AppendedPotentialDerivative(derivativeList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PotentialDerivativeInterface.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AppendedPotentialDerivative.class;
    }
}
