package dr.inferencexml.operators;

import dr.inference.distribution.IndependentNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledParameter;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.inferencexml.distribution.IndependentNormalDistributionModelParser;
import dr.inferencexml.model.ProductParameterParser;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class MultiplicativeGammaGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String MULTIPLICATIVE_PROVIDER = "multiplicativeGammaGibbsProvider";

    private static final String PRECISION = "precisionType";
    private static final String LOCAL = "local";
    private static final String GLOBAL = "global";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        IndependentNormalDistributionModel normalDistribution = (IndependentNormalDistributionModel) xo.getChild(IndependentNormalDistributionModel.class);

        if (!normalDistribution.precisionUsed()) {
            throw new XMLParseException("The " +
                    IndependentNormalDistributionModelParser.INDEPENDENT_NORMAL_DISTRIBUTION_MODEL +
                    " object must be constructed with a precision, not variance.");
        }

        Parameter data = normalDistribution.getData();
        Parameter mean = normalDistribution.getMean();

        Parameter precision = normalDistribution.getPrecision();

        if (!(precision instanceof ScaledParameter)) {
            throw new XMLParseException("The precision parameter in the " +
                    IndependentNormalDistributionModelParser.INDEPENDENT_NORMAL_DISTRIBUTION_MODEL +
                    " object must be a " + ProductParameterParser.PRODUCT_PARAMETER + " with a '" +
                    ProductParameterParser.SCALE + "' component.");
        }

        ScaledParameter scaledPrecision = (ScaledParameter) precision;

        Parameter globalPrecision = scaledPrecision.getScaleParam();
        Parameter localPrecision = scaledPrecision.getVecParam();

        String precisionType = xo.getStringAttribute(PRECISION);
        if (precisionType.equalsIgnoreCase(LOCAL)) {

            return new GammaGibbsProvider.LocalMultiplicativeGammaGibbsProvider(data, mean,
                    globalPrecision, localPrecision);

        } else if (precisionType.equalsIgnoreCase(GLOBAL)) {

            return new GammaGibbsProvider.GlobalMultiplicativeGammaGibbsProvider(data, mean,
                    globalPrecision, localPrecision);

        } else {
            throw new XMLParseException("The '" + PRECISION + "' attribute must be either '" + LOCAL + "' or '" +
                    GLOBAL + "'.");
        }


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule(PRECISION),
                new ElementRule(IndependentNormalDistributionModel.class)
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
