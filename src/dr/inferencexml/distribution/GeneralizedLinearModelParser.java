package dr.inferencexml.distribution;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LinearRegression;
import dr.inference.distribution.LogLinearModel;
import dr.inference.distribution.LogisticRegression;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class GeneralizedLinearModelParser extends AbstractXMLObjectParser {

    public static final String GLM_LIKELIHOOD = "glmModel";

    public static final String DEPENDENT_VARIABLES = "dependentVariables";
    public static final String INDEPENDENT_VARIABLES = "independentVariables";
    public static final String BASIS_MATRIX = "basis";
    public static final String FAMILY = "family";
    public static final String SCALE_VARIABLES = "scaleVariables";
    public static final String INDICATOR = "indicator";
    public static final String LOGISTIC_REGRESSION = "logistic";
    public static final String NORMAL_REGRESSION = "normal";
    public static final String LOG_NORMAL_REGRESSION = "logNormal";
    public static final String LOG_LINEAR = "logLinear";
//    public static final String LOG_TRANSFORM = "logDependentTransform";
    public static final String RANDOM_EFFECTS = "randomEffects";
    public static final String CHECK_IDENTIFIABILITY = "checkIdentifiability";

    public String getParserName() {
        return GLM_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DEPENDENT_VARIABLES);
        Parameter dependentParam = null;
        if (cxo != null)
            dependentParam = (Parameter) cxo.getChild(Parameter.class);

        String family = xo.getStringAttribute(FAMILY);
        GeneralizedLinearModel glm;
        if (family.compareTo(LOGISTIC_REGRESSION) == 0) {
            glm = new LogisticRegression(dependentParam);
        } else if (family.compareTo(NORMAL_REGRESSION) == 0) {
            glm = new LinearRegression(dependentParam, false);
        } else if (family.compareTo(LOG_NORMAL_REGRESSION) == 0) {
            glm = new LinearRegression(dependentParam, true);
        } else if (family.compareTo(LOG_LINEAR) == 0) {
            glm = new LogLinearModel(dependentParam);
        } else
            throw new XMLParseException("Family '" + family + "' is not currently implemented");

        if (glm.requiresScale()) {
            cxo = xo.getChild(SCALE_VARIABLES);
            Parameter scaleParameter = null;
//                DesignMatrix designMatrix = null;
            Parameter scaleDesign = null;
            if (cxo != null) {
                scaleParameter = (Parameter) cxo.getChild(Parameter.class);
                XMLObject gxo = cxo.getChild(INDICATOR);
                if (gxo != null)
                    scaleDesign = (Parameter) gxo.getChild(Parameter.class);
//                    designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
            }
            if (scaleParameter == null)
                throw new XMLParseException("Family '" + family + "' requires scale parameters");
            if (scaleDesign == null)
                scaleDesign = new Parameter.Default(dependentParam.getDimension(), 0.0);
            else {
                if (scaleDesign.getDimension() != dependentParam.getDimension())
                    throw new XMLParseException("Scale and scaleDesign parameters must be the same dimension");
                for (int i = 0; i < scaleDesign.getDimension(); i++) {
                    double value = scaleDesign.getParameterValue(i);
                    if (value < 1 || value > scaleParameter.getDimension())
                        throw new XMLParseException("Invalid scaleDesign value");
                    scaleDesign.setParameterValue(i, value - 1);
                }
            }

            glm.addScaleParameter(scaleParameter, scaleDesign);
        }

        addIndependentParameters(xo, glm, dependentParam);
        addRandomEffects(xo, glm, dependentParam);

        boolean checkIdentifiability = xo.getAttribute(CHECK_IDENTIFIABILITY, true);
        if (checkIdentifiability) {
            if (!glm.getAllIndependentVariablesIdentifiable()) {
                throw new XMLParseException("All design matrix predictors are not identifiable in "+  xo.getId());
            }
        }

        return glm;
    }

    public void addRandomEffects(XMLObject xo, GeneralizedLinearModel glm,
                                 Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(RANDOM_EFFECTS) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter randomEffect = (Parameter) cxo.getChild(Parameter.class);
                checkRandomEffectsDimensions(randomEffect, dependentParam);
                glm.addRandomEffectsParameter(randomEffect);
            }
        }
    }

    public void addIndependentParameters(XMLObject xo, GeneralizedLinearModel glm,
                                         Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(INDEPENDENT_VARIABLES) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
                DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
                checkDimensions(independentParam, dependentParam, designMatrix);
                cxo = cxo.getChild(INDICATOR);
                Parameter indicator = null;
                if (cxo != null) {
                    indicator = (Parameter) cxo.getChild(Parameter.class);
                    if (indicator.getDimension() != independentParam.getDimension())
                        throw new XMLParseException("dim(" + independentParam.getId() + ") != dim(" + indicator.getId() + ")");
                }
                checkFullRank(designMatrix);
                glm.addIndependentParameter(independentParam, designMatrix, indicator);
            }
        }
    }

    private void checkFullRank(DesignMatrix designMatrix) throws XMLParseException {
        int fullRank = designMatrix.getColumnDimension();

        SingularValueDecomposition svd = new SingularValueDecomposition(
                new DenseDoubleMatrix2D(designMatrix.getParameterAsMatrix()));
        int realRank = svd.rank();
        if (realRank != fullRank) {
            throw new XMLParseException(
                "rank(" + designMatrix.getId() + ") = " + realRank +
                        ".\nMatrix is not of full rank as colDim(" + designMatrix.getId() + ") = " + fullRank        
            );
        }
    }

    private void checkRandomEffectsDimensions(Parameter randomEffect, Parameter dependentParam)
            throws XMLParseException {
        if (dependentParam != null) {
            if (randomEffect.getDimension() != dependentParam.getDimension()) {
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + randomEffect.getId() + ")"
                );
            }
        }
    }

    private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
            throws XMLParseException {
        if (dependentParam != null) {
            if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                    (independentParam.getDimension() != designMatrix.getColumnDimension()))
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
                );
        } else {
            if (independentParam.getDimension() != designMatrix.getColumnDimension()) {
                throw new XMLParseException(
                        "dim(" + independentParam.getId() + ") is incompatible with dim (" + designMatrix.getId() + ")"
                );
            }
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FAMILY),
            AttributeRule.newBooleanRule(CHECK_IDENTIFIABILITY, true),
            new ElementRule(DEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(INDEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true),
                            new ElementRule(DesignMatrix.class),
                            new ElementRule(INDICATOR,
                                    new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }, true),
                    }, 1, 3),
            new ElementRule(RANDOM_EFFECTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, 0, 3),
//				new ElementRule(BASIS_MATRIX,
//						new XMLSyntaxRule[]{new ElementRule(DesignMatrix.class)})
    };

    public String getParserDescription() {
        return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }
}
