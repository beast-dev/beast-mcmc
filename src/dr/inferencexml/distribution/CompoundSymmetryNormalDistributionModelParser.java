package dr.inferencexml.distribution;

import dr.inference.distribution.CompoundSymmetryNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class CompoundSymmetryNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "compoundSymmetryNormalDistributionModel";
    private static final String DIMENSION = "dim";
    private static final String MARGINAL_VARIANCE = "variance";
    private static final String CORRELATION = "rho";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        XMLObject cxo = xo.getChild(MARGINAL_VARIANCE);
        Parameter variance = (Parameter) cxo.getChild(Parameter.class);

        if (variance.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("variance must be > 0.0");
        }

        cxo = xo.getChild(CORRELATION);
        Parameter rho = (Parameter) cxo.getChild(Parameter.class);

        if (Math.abs(rho.getParameterValue(0)) >= 1.0) {
            throw new XMLParseException("|Rho| must be < 1.0");
        }
        return new CompoundSymmetryNormalDistributionModel(dim, variance, rho);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(MARGINAL_VARIANCE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(CORRELATION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundSymmetryNormalDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }
}
