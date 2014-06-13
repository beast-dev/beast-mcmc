package dr.inferencexml.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.DiagonalMatrix;
import dr.math.Polynomial;
import dr.inference.operators.FactorGibbsOperator;
import dr.xml.*;

/**
 @author Max Tolkoff
 */
public class FactorGibbsOperatorParser extends AbstractXMLObjectParser {
    private final String FACTOR_GIBBS_SAMPLER="factorGibbsOperator";
    private final String WEIGHT="weight";
    private final String RANDOM_SCAN="randomScan";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp= (String) xo.getAttribute(WEIGHT);
        double weight=Double.parseDouble(weightTemp);
        DiagonalMatrix diffusionMatrix;
        diffusionMatrix=(DiagonalMatrix) xo.getChild(DiagonalMatrix.class);
        LatentFactorModel LFM =(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        boolean randomScan=xo.getAttribute(RANDOM_SCAN, true);
        return new FactorGibbsOperator(LFM, weight, randomScan, diffusionMatrix);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
//            new ElementRule(CompoundParameter.class),
            new ElementRule(DiagonalMatrix.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FactorGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_GIBBS_SAMPLER;
    }
}
