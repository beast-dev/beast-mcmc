package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.MultivariateOUModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.model.DesignMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.GeneralizedLinearModelParser;
import dr.xml.*;

/**
 *
 */
public class MultivariateOUModelParser extends AbstractXMLObjectParser {

    public static final String MVOU_MODEL = "multivariateOUModel";
    public static final String MVOU_TYPE = "MVOU";
    public static final String DATA = "data";
    public static final String TIME = "times";
    public static final String DESIGN = "design";

    public String getParserName() {
        return MVOU_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        Parameter effectParameter = (Parameter) ((XMLObject) xo.getChild(DATA)).getChild(Parameter.class);
        Parameter timesParameter = (Parameter) ((XMLObject) xo.getChild(TIME)).getChild(Parameter.class);
        Parameter designParameter = (Parameter) ((XMLObject) xo.getChild(DESIGN)).getChild(Parameter.class);
        MatrixParameter gammaParameter = (MatrixParameter) xo.getChild(MatrixParameter.class);

        if (effectParameter.getDimension() != timesParameter.getDimension() ||
                effectParameter.getDimension() != designParameter.getDimension()) {
//				System.err.println("dim(effect) " +effectParameter.getDimension());
//				System.err.println("dim(times) "+timesParameter.getDimension());
//				System.err.println("dim(design) "+designParameter.getDimension());
            throw new XMLParseException("dim(" + effectParameter.getStatisticName() +
                    ") != dim(" + timesParameter.getStatisticName() + ") != dim(" + designParameter.getStatisticName() +
                    ") in " + xo.getName() + " element");
        }

        MultivariateOUModel glm = new MultivariateOUModel(substitutionModel, effectParameter, gammaParameter,
                timesParameter.getParameterValues(), designParameter.getParameterValues());

        addIndependentParameters(xo, glm, effectParameter);

        // todo Confirm that design vector is consistent with substitution model
        // todo Confirm that design vector is ordered 1,\ldots,K,1,\ldots,K, etc.

        return glm;

    }

    public void addIndependentParameters(XMLObject xo, GeneralizedLinearModel glm,
                                         Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(GeneralizedLinearModelParser.INDEPENDENT_VARIABLES) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
                DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
                checkDimensions(independentParam, dependentParam, designMatrix);
                glm.addIndependentParameter(independentParam, designMatrix, null);
            }
        }
    }

    private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
            throws XMLParseException {
        if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                (independentParam.getDimension() != designMatrix.getColumnDimension())) {
            System.err.println(dependentParam.getDimension());
            System.err.println(independentParam.getDimension());
            System.err.println(designMatrix.getRowDimension() + " rows");
            System.err.println(designMatrix.getColumnDimension() + " cols");
            throw new XMLParseException(
                    "dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
            );
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SubstitutionModel.class),
            new ElementRule(MatrixParameter.class),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)}),
            new ElementRule(TIME, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)}),
            new ElementRule(DESIGN, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)}),
            new ElementRule(GeneralizedLinearModelParser.INDEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, 0, 3),
    };

    public String getParserDescription() {
        return "Describes a multivariate OU process";
    }

    public Class getReturnType() {
        return MultivariateOUModel.class;
    }
}
