package dr.inference.distribution;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.MultivariateFunction;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 */
public abstract class GeneralizedLinearModel extends AbstractModel implements Likelihood, MultivariateFunction {
//		, RealFunctionOfSeveralVariablesWithGradient {

    public static final String GLM_LIKELIHOOD = "glmModel";

    public static final String DEPENDENT_VARIABLES = "dependentVariables";
    public static final String INDEPENDENT_VARIABLES = "independentVariables";
    public static final String BASIS_MATRIX = "basis";
    public static final String FAMILY = "family";
    public static final String SCALE_VARIABLES = "scaleVariables";
    public static final String LOGISTIC_REGRESSION = "logistic";
    public static final String NORMAL_REGRESSION = "normal";
    public static final String LOG_LINEAR = "logLinear";
//	public static final String RANDOM_EFFECTS = "randomEffects";

    protected Parameter dependentParam;
    protected List<Parameter> independentParam;
    protected List<double[][]> designMatrix; // fixed constants, access as double[][] to save overhead

    protected double[][] scaleDesignMatrix;
    protected Parameter scaleParameter;

    protected int numIndependentVariables = 0;
    protected int N;

    protected List<Parameter> randomEffects;

    public GeneralizedLinearModel(Parameter dependentParam) { //, Parameter independentParam,
//	                              DesignMatrix designMatrix) {
        super(GLM_LIKELIHOOD);
        this.dependentParam = dependentParam;
//		this.independentParam = independentParam;
//		this.designMatrix = designMatrix;
//		addParameter(independentParam);
        if (dependentParam != null) {
            addParameter(dependentParam);
            N = dependentParam.getDimension();
        } else
            N = 0;
    }

    public void addIndependentParameter(Parameter effect, DesignMatrix matrix) {
        if (designMatrix == null)
            designMatrix = new ArrayList<double[][]>();
        if (independentParam == null)
            independentParam = new ArrayList<Parameter>();
        if (N == 0) {
            N = matrix.getRowDimension();
        }
        designMatrix.add(matrix.getParameterAsMatrix());
        independentParam.add(effect);
        if (designMatrix.size() != independentParam.size())
            throw new RuntimeException("Independent variables and their design matrices are out of sync");
        addParameter(effect);
        numIndependentVariables++;
        System.out.println("\tAdding independent predictors '" + effect.getStatisticName() + "' with design matrix '" + matrix.getStatisticName() + "'");
    }

    public int getNumberOfEffects() {
        return numIndependentVariables;
    }

    public double[] getXBeta() {

        double[] xBeta = new double[N];

        for (int j = 0; j < numIndependentVariables; j++) {
            Parameter beta = independentParam.get(j);
            double[][] X = designMatrix.get(j);
            final int K = beta.getDimension();
            for (int k = 0; k < K; k++) {
                final double betaK = beta.getParameterValue(k);
                for (int i = 0; i < N; i++)
                    xBeta[i] += X[i][k] * betaK;
            }
        }

        return xBeta;

    }

    public Parameter getEffect(int j) {
        return independentParam.get(j);
    }

    public Parameter getDependentVariable() {
        return dependentParam;
    }

    public double[] getXBeta(int j) {

        double[] xBeta = new double[N];

        Parameter beta = independentParam.get(j);
        double[][] X = designMatrix.get(j);
        final int K = beta.getDimension();
        for (int k = 0; k < K; k++) {
            final double betaK = beta.getParameterValue(k);
            for (int i = 0; i < N; i++)
                xBeta[i] += X[i][k] * betaK;
        }

        return xBeta;

    }

    public int getEffectNumber(Parameter effect) {
        return independentParam.indexOf(effect);
    }

//	public double[][] getXtScaleX(int j) {
//
//		final Parameter beta = independentParam.get(j);
//		double[][] X = designMatrix.get(j);
//		final int dim = X[0].length;
//
//		if( dim != beta.getDimension() )
//			throw new RuntimeException("should have checked eariler");
//
//		double[] scale = getScale();
//
//
//	}

    public double[][] getX(int j) {
        return designMatrix.get(j);
    }


    public double[] getScale() {

        double[] scale = new double[N];

        final int K = scaleParameter.getDimension();
        for (int k = 0; k < K; k++) {
            final double scaleK = scaleParameter.getParameterValue(k);
            for (int i = 0; i < N; i++)
                scale[i] += scaleDesignMatrix[i][k] * scaleK;
        }

        return scale;
    }


    public double[][] getScaleAsMatrix() {

        double[][] scale = new double[N][N];

        return scale;
    }

//	protected abstract double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient);

    protected abstract double calculateLogLikelihood(double[] beta);

    protected abstract double calculateLogLikelihood();

    protected abstract boolean confirmIndependentParameters();

    protected abstract boolean requiresScale();

    private void addScaleParameter(Parameter scaleParameter, DesignMatrix matrix) {
        this.scaleParameter = scaleParameter;
        this.scaleDesignMatrix = matrix.getParameterAsMatrix();
        addParameter(scaleParameter);
    }

/*	// **************************************************************
	// RealFunctionOfSeveralVariablesWithGradient IMPLEMENTATION
	// **************************************************************


	public double eval(double[] beta, double[] gradient) {
		return calculateLogLikelihoodAndGradient(beta, gradient);
	}


	public double eval(double[] beta) {
		return calculateLogLikelihood(beta);
	}


	public int getNumberOfVariables() {
		return independentParam.getDimension();
	}*/

    // ************
    //       Mutlivariate implementation
    // ************


    public double evaluate(double[] beta) {
        return calculateLogLikelihood(beta);
    }

    public int getNumArguments() {
        int total = 0;
        for (Parameter effect : independentParam)
            total += effect.getDimension();
        return total;
    }

    public double getLowerBound(int n) {
        int which = n;
        int k = 0;
        while (which > independentParam.get(k).getDimension()) {
            which -= independentParam.get(k).getDimension();
            k++;
        }
        return independentParam.get(k).getBounds().getLowerLimit(which);
    }

    public double getUpperBound(int n) {
        int which = n;
        int k = 0;
        while (which > independentParam.get(k).getDimension()) {
            which -= independentParam.get(k).getDimension();
            k++;
        }
        return independentParam.get(k).getBounds().getUpperLimit(which);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
//        fireModelChanged();
    }

    protected void storeState() {
        // No internal states to save
    }

    protected void restoreState() {
        // No internal states to restore
    }

    protected void acceptState() {
        // Nothing to do
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }


    @Override
    public String toString() {
        return super.toString() + ": " + getLogLikelihood();
    }

    public void makeDirty() {
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

//    /**
//     * @return the log columns.
//     */
//    public LogColumn[] getColumns() {
//        return new dr.inference.loggers.LogColumn[]{
//                new LikelihoodColumn(getId())
//        };
//    }
//
//    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
//        public LikelihoodColumn(String label) {
//            super(label);
//        }
//
//        public double getDoubleValue() {
//            return getLogLikelihood();
//        }
//    }

    public LogColumn[] getColumns() {
        LogColumn[] output = new LogColumn[N];
        for(int i=0; i<N; i++)
            output[i] = new NumberArrayColumn(getId()+i,i);
        return output;
    }

    private class NumberArrayColumn extends NumberColumn {

        private int index;

        public NumberArrayColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        public double getDoubleValue() {
            return getXBeta()[index];
        }
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }


    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GLM_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(DEPENDENT_VARIABLES);
            Parameter dependentParam = null;
            if (cxo != null)
                dependentParam = (Parameter) cxo.getChild(Parameter.class);

            String family = xo.getStringAttribute(FAMILY);
            GeneralizedLinearModel glm;
            if (family.compareTo(LOGISTIC_REGRESSION) == 0) {
                glm = new LogisticRegression(dependentParam);
            } else if (family.compareTo(NORMAL_REGRESSION) == 0) {
                glm = new LinearRegression(dependentParam);
            } else if (family.compareTo(LOG_LINEAR) == 0) {
                glm = new LogLinearModel(dependentParam);
            } else
                throw new XMLParseException("Family '" + family + "' is not currently implemented");

            if (glm.requiresScale()) {
                cxo = (XMLObject) xo.getChild(SCALE_VARIABLES);
                Parameter scaleParameter = null;
                DesignMatrix designMatrix = null;
                if (cxo != null) {
                    scaleParameter = (Parameter) cxo.getChild(Parameter.class);
                    designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
                }
                if (scaleParameter == null || designMatrix == null)
                    throw new XMLParseException("Family '" + family + "' requires scale parameters");
                checkDimensions(scaleParameter, dependentParam, designMatrix);
                glm.addScaleParameter(scaleParameter, designMatrix);
            }

            addIndependentParameters(xo, glm, dependentParam);

            return glm;
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
                    glm.addIndependentParameter(independentParam, designMatrix);
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
                           "dim(" +independentParam.getId() + ") is incompatible with dim (" + designMatrix.getId() +")"
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(FAMILY),
                new ElementRule(DEPENDENT_VARIABLES,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
                new ElementRule(INDEPENDENT_VARIABLES,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, 1, 3),
//				new ElementRule(BASIS_MATRIX,
//						new XMLSyntaxRule[]{new ElementRule(DesignMatrix.class)})
        };

        public String getParserDescription() {
            return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };


}
