package dr.inferencexml.model;

import dr.inference.model.CompoundGradientProvider;
import dr.inference.model.GradientProvider;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivativeParser extends AbstractXMLObjectParser {

    public final static String SUM_DERIVATIVE = "appendedPotentialDerivative";
    public static final String SUM_DERIVATIVE2 = "compoundGradientForLikelihood";

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

        List<GradientProvider> derivativeList = new ArrayList<GradientProvider>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            derivativeList.add((GradientProvider) xo.getChild(i));
        }


        return new CompoundGradientProvider(derivativeList);
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
        return CompoundGradientProvider.class;
    }
}
