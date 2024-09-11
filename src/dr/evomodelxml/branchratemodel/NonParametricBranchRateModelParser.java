package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.NonParametricBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.xml.*;



public class NonParametricBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "nonParametricRates";
    private static final String DEGREE = "degree";
    private static final String COEFFICIENTS = "coefficients";
    private static final String BOUNDARYFACTOR = "boundaryFactor";
    private static final String MARGINALVARIANCE = "marginalVariance";
    private static final String LENGTHSCALE = "lengthScale";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        double degree = xo.getDoubleAttribute(DEGREE);
        Parameter coefficients = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        Parameter boundaryFactor = (Parameter) xo.getElementFirstChild(BOUNDARYFACTOR);
        Parameter marginalVariance = (Parameter) xo.getElementFirstChild(MARGINALVARIANCE);
        Parameter lengthScale = (Parameter) xo.getElementFirstChild(LENGTHSCALE);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new NonParametricBranchRateModel(id, tree, degree, coefficients, boundaryFactor, marginalVariance, lengthScale);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NonParametricBranchRateModelParser.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            AttributeRule.newDoubleRule(DEGREE),
            new ElementRule(COEFFICIENTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(BOUNDARYFACTOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(MARGINALVARIANCE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(LENGTHSCALE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
    };
}
