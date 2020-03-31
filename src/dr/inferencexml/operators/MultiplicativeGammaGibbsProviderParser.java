package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class MultiplicativeGammaGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String MULTIPLICATIVE_PROVIDER = "multiplicativeGammaGibbsProvider";
    private static final String MEAN = "mean";
    private static final String GLOBAL_PREC = "globalPrecision";
    private static final String LOCAL_PREC = "localPrecision";

    private static final String PRECISION = "precisionType";
    private static final String LOCAL = "local";
    private static final String GLOBAL = "global";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        XMLObject mxo = xo.getChild(MEAN);
        Parameter mean = (Parameter) mxo.getChild(Parameter.class);

        XMLObject gxo = xo.getChild(GLOBAL_PREC);
        Parameter globalPrecision = (Parameter) gxo.getChild(Parameter.class);

        XMLObject lxo = xo.getChild(LOCAL_PREC);
        Parameter localPrecision = (Parameter) lxo.getChild(Parameter.class);

        String precisionType = xo.getStringAttribute(PRECISION);
        if (precisionType.equalsIgnoreCase(LOCAL)) {
            return new GammaGibbsProvider.LocalMultiplicativeGammaGibbsProvider(parameter, mean,
                    globalPrecision, localPrecision);
        } else if (precisionType.equalsIgnoreCase(GLOBAL)) {
            return new GammaGibbsProvider.GlobalMultiplicativeGammaGibbsProvider(parameter, mean,
                    globalPrecision, localPrecision);
        } else {
            throw new XMLParseException("Attribute '" + PRECISION + "' must be either '" + GLOBAL + "' or '" +
                    LOCAL + "'.");
        }


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule(PRECISION),
                new ElementRule(Parameter.class),
                new ElementRule(MEAN,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }),

                new ElementRule(GLOBAL_PREC,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }),

                new ElementRule(LOCAL_PREC,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        })


        };
    }

    @Override
    public String getParserDescription() {
        return "Returns sufficient statistics for a normally distributed parameter with" +
                " multiplicative gamma precision prior";
    }

    @Override
    public Class getReturnType() {
        return GammaGibbsProvider.class;
    }

    @Override
    public String getParserName() {
        return MULTIPLICATIVE_PROVIDER;
    }
}
