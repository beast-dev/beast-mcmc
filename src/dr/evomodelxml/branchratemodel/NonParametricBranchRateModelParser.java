package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.NonParametricBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.xml.*;



public class NonParametricBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "nonParametricRates";
    private static final String LASTSAMPLINGTIME = "lastSamplingTime";
    private static final String COEFFICIENTS = "coefficients";
    private static final String BOUNDARY = "boundary";
    private static final String SCALEFACTOR = "scaleFactor";
    private static final String LENGTHSCALE = "lengthScale";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter lastSamplingTime = (Parameter) xo.getElementFirstChild(LASTSAMPLINGTIME);
        Parameter coefficients = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        Parameter boundary = (Parameter) xo.getElementFirstChild(BOUNDARY);
        Parameter scaleFactor = (Parameter) xo.getElementFirstChild(SCALEFACTOR);
        Parameter lengthScale = (Parameter) xo.getElementFirstChild(LENGTHSCALE);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new NonParametricBranchRateModel(id, tree, lastSamplingTime, coefficients, boundary, scaleFactor, lengthScale);
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
            new ElementRule(LASTSAMPLINGTIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(COEFFICIENTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(BOUNDARY,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(SCALEFACTOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(LENGTHSCALE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
    };
}
