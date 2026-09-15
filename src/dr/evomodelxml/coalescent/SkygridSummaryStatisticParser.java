package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFGaussianApproximation;
import dr.evomodel.coalescent.SkygridSummaryStatistic;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

public class SkygridSummaryStatisticParser extends AbstractXMLObjectParser {

    public static final String SKYGRID_SUMMARY_STATISTIC = "skygridSummaryStatistic";
    public static final String TYPE = "type";
    public static final String COEFFICIENT = "coefficient";
    public static final String BETA_PRIOR_PRECISION = "betaPriorPrecision";
    public static final String BETA_PRIOR_MEAN = "betaPriorMean";
    public static final String TAU_SHAPE = "tauShape";
    public static final String TAU_RATE = "tauRate";

    public String getParserName(){
        return SKYGRID_SUMMARY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GMRFGaussianApproximation approximation =
                (GMRFGaussianApproximation) xo.getChild(GMRFGaussianApproximation.class);

        if(approximation == null){
            throw new XMLParseException("SkygridSummaryStatistic requires a GMRFGaussianApproximation");
        }

        String typeString = xo.getStringAttribute(TYPE);
        SkygridSummaryStatistic.Type type = parseType(typeString);

        Integer coefficientIndex = null;
        if(xo.hasAttribute(COEFFICIENT)){
            coefficientIndex = xo.getIntegerAttribute(COEFFICIENT)-1;
            if(coefficientIndex < 0){
                throw new XMLParseException("invalid coefficient index: " + coefficientIndex);
            }
        }

        double[][] betaPriorPrecision = null;
        if(xo.hasChildNamed(BETA_PRIOR_PRECISION)){
            XMLObject cxo = xo.getChild(BETA_PRIOR_PRECISION);
            MatrixParameter matrixParameter = (MatrixParameter) cxo.getChild(MatrixParameter.class);
            betaPriorPrecision = convertToDoubleArray(matrixParameter);
        }

        double[] betaPriorMean = null;
        if(xo.hasChildNamed(BETA_PRIOR_MEAN)){
            XMLObject cxo = xo.getChild(BETA_PRIOR_MEAN);
            Parameter meanParameter = (Parameter) cxo.getChild(Parameter.class);
            betaPriorMean = meanParameter.getParameterValues();
        }

        Double tauShape = null;
        if(xo.hasChildNamed(TAU_SHAPE)){
            XMLObject cxo = xo.getChild(TAU_SHAPE);
            Parameter tauShapeParameter = (Parameter) cxo.getChild(Parameter.class);
            tauShape = tauShapeParameter.getParameterValue(0);
        }

        Double tauRate = null;
        if(xo.hasChildNamed(TAU_RATE)){
            XMLObject cxo = xo.getChild(TAU_RATE);
            Parameter tauRateParameter = (Parameter) cxo.getChild(Parameter.class);
            tauRate = tauRateParameter.getParameterValue(0);
        }

        switch(type){
            case MUTUAL_INFORMATION:
            case BETA_GIVEN_GAMMA_VARIANCE_DIAG:
            case BETA_TILDE_VARIANCE_DIAG:
                requireXML(betaPriorPrecision != null, typeString, "a betaPriorPrecision child element");
                break;
            case BETA_GIVEN_GAMMA_MEAN:
                requireXML(betaPriorPrecision != null, typeString, "a betaPriorPrecision child element");
                requireXML(betaPriorMean != null, typeString, "a betaPriorMean chld element");
                break;
            case BETA_TILDE_MEAN:
                requireXML(betaPriorPrecision != null, typeString, "a betaPriorPrecision child element");
                requireXML(betaPriorMean != null, typeString, "a betaPriorMean chld element");
                break;
                case TAU_CONDITIONAL_SHAPE:
                requireXML(tauShape != null, typeString, "a tauShape child element");
                break;
            case TAU_CONDITIONAL_RATE:
                requireXML(tauRate != null, typeString, "a tauRate child element");
                break;
                default:
                    break;
        }

        Parameter gamma = null;
        switch(type){
            case LAGRANGE_BOUND:
            case BETA_GIVEN_GAMMA_MEAN:
            case TAU_CONDITIONAL_RATE:
                gamma = approximation.getLikelihood().getPopSizeParameter();
                break;
                default:
                    break;
        }

        String id = xo.hasId() ? xo.getId() : typeString;
        return new SkygridSummaryStatistic(id, approximation, type, coefficientIndex, betaPriorPrecision,
                betaPriorMean, tauShape, tauRate, gamma);

    }

    private void requireXML(boolean condition, String typeString, String message) throws XMLParseException {
        if(!condition){
            throw new XMLParseException("type=\"" + typeString + "\" requires " + message + ".");
        }
    }

    private SkygridSummaryStatistic.Type parseType(String typeString) throws XMLParseException {
        switch (typeString) {
            case "retainedInformationRatio":
                return SkygridSummaryStatistic.Type.RETAINED_INFORMATION_RATIO;
            case "mutualInformation":
                return SkygridSummaryStatistic.Type.MUTUAL_INFORMATION;
            case "lagrangeBound":
                return SkygridSummaryStatistic.Type.LAGRANGE_BOUND;
            case "gammaHat":
                return SkygridSummaryStatistic.Type.GAMMA_HAT;
            case "sigmaGInverseDiag":
                return SkygridSummaryStatistic.Type.SIGMA_G_INVERSE_DIAG;
            case "gammaTildeMean":
                return SkygridSummaryStatistic.Type.GAMMA_TILDE_MEAN;
            case "gammaTildeVarianceDiag":
                return SkygridSummaryStatistic.Type.GAMMA_TILDE_VARIANCE_DIAG;
            case "betaGivenGammaMean":
                return SkygridSummaryStatistic.Type.BETA_GIVEN_GAMMA_MEAN;
            case "betaGivenGammaVarianceDiag":
                return SkygridSummaryStatistic.Type.BETA_GIVEN_GAMMA_VARIANCE_DIAG;
            case "betaTildeMean":
                return SkygridSummaryStatistic.Type.BETA_TILDE_MEAN;
            case "betaTildeVarianceDiag":
                return SkygridSummaryStatistic.Type.BETA_TILDE_VARIANCE_DIAG;
            case "tauConditionalShape":
                return SkygridSummaryStatistic.Type.TAU_CONDITIONAL_SHAPE;
            case "tauConditionalRate":
                return SkygridSummaryStatistic.Type.TAU_CONDITIONAL_RATE;
            default:
                throw new XMLParseException("Unrecognized skygridSummaryStatistic: " + typeString);
        }
    }

    private double[][] convertToDoubleArray(MatrixParameter matrixParameter){
        int rows = matrixParameter.getRowDimension();
        int columns = matrixParameter.getColumnDimension();
        double[][] result = new double[rows][columns];
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                result[i][j] = matrixParameter.getParameterValue(i, j);
            }
        }
        return result;
    }

    public String getParserDescription(){
        return "Computes summary statistics for skygrid-GLM analysis.";
    }

    public Class getReturnType(){
        return SkygridSummaryStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules(){
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TYPE),
            AttributeRule.newIntegerRule(COEFFICIENT, true),
            new ElementRule(GMRFGaussianApproximation.class),
            new ElementRule(BETA_PRIOR_PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }, true),
            new ElementRule(BETA_PRIOR_MEAN, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(TAU_SHAPE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(TAU_RATE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };
}