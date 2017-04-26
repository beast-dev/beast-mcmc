package dr.inferencexml.model;

import dr.inference.model.GradientProvider;
import dr.inference.model.SumDerivative;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

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
        List<GradientProvider> derivativeList = new ArrayList<GradientProvider>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            derivativeList.add((GradientProvider) xo.getChild(i));
        }


        return new SumDerivative(derivativeList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientProvider.class, 1, Integer.MAX_VALUE),
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
