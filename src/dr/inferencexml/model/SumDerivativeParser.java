package dr.inferencexml.model;

import dr.inference.model.GradientWrtParameterProvider;
import dr.inference.model.SumDerivative;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class SumDerivativeParser extends AbstractXMLObjectParser{
    public final static String SUM_DERIVATIVE = "sumDerivative";
    public final static String SUM_DERIVATIVE2 = "jointGradient";


    @Override
    public String getParserName() {
        return SUM_DERIVATIVE;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { SUM_DERIVATIVE, SUM_DERIVATIVE2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<GradientWrtParameterProvider> derivativeList = new ArrayList<GradientWrtParameterProvider>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            GradientWrtParameterProvider grad = (GradientWrtParameterProvider) xo.getChild(i);
            derivativeList.add(grad);
        }

        return new SumDerivative(derivativeList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class, 1, Integer.MAX_VALUE),
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
