package dr.inferencexml.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.model.DiagonalMatrix;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.FactorOperator;
import dr.xml.*;

/**
 @author Max Tolkoff
 */
public class FactorOperatorParser extends AbstractXMLObjectParser {
    private final String FACTOR_OPERATOR="factorOperator";
    private final String WEIGHT="weight";
    private final String RANDOM_SCAN="randomScan";
    private final String SCALE_FACTOR="scaleFactor";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        CoercionMode mode = CoercionMode.parseMode(xo);
        String scaleFactorTemp=(String)xo.getAttribute(SCALE_FACTOR);
        double scaleFactor=Double.parseDouble(scaleFactorTemp);
        String weightTemp= (String) xo.getAttribute(WEIGHT);
        double weight=Double.parseDouble(weightTemp);
        DiagonalMatrix diffusionMatrix;
        diffusionMatrix=(DiagonalMatrix) xo.getChild(DiagonalMatrix.class);
        LatentFactorModel LFM =(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        boolean randomScan=xo.getAttribute(RANDOM_SCAN, true);
        return new FactorOperator(LFM, weight, randomScan, diffusionMatrix, scaleFactor, mode);
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
            AttributeRule.newDoubleRule(SCALE_FACTOR),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FactorOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_OPERATOR;
    }
}
