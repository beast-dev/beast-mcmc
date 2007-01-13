package dr.inference.distribution;

import dr.inference.loggers.LogColumn;
import dr.inference.model.*;
import dr.math.MultivariateFunction;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Dec 29, 2006
 * Time: 11:01:02 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GeneralizedLinearModel extends AbstractModel implements Likelihood, MultivariateFunction {
//		, RealFunctionOfSeveralVariablesWithGradient {

	public static final String GLM_LIKELIHOOD = "glmLikelihood";

	public static final String DEPENDENT_VARIABLES = "dependentVariables";
	public static final String INDEPENDENT_VARIABLES = "independentVariables";
	public static final String BASIS_MATRIX = "basis";
	public static final String FAMILY = "family";
	public static final String SCALE = "scaleParameter";
	public static final String LOGISTIC_REGRESSION = "logistic";

	protected Parameter dependentParam;
	protected Parameter independentParam;
	protected DesignMatrix designMatrix;
	protected Parameter scaleParameter;

	public GeneralizedLinearModel(Parameter dependentParam, Parameter independentParam,
	                              DesignMatrix designMatrix) {
		super(GLM_LIKELIHOOD);
		this.dependentParam = dependentParam;
		this.independentParam = independentParam;
		this.designMatrix = designMatrix;
		addParameter(independentParam);
		addParameter(dependentParam);
	}

//	protected abstract double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient);

	protected abstract double calculateLogLikelihood(double[] beta);

	protected abstract double calculateLogLikelihood();

	protected abstract boolean confirmIndependentParameters();

	protected abstract boolean requiresScale();

	private void addScaleParameter(Parameter scaleParameter) {
		this.scaleParameter = scaleParameter;
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
		return independentParam.getDimension();
	}

	public double getLowerBound(int n) {
		return independentParam.getBounds().getLowerLimit(n);
	}

	public double getUpperBound(int n) {
		return independentParam.getBounds().getUpperLimit(n);
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {

	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
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
		return calculateLogLikelihood(independentParam.getParameterValues());
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

	/**
	 * @return the log columns.
	 */
	public LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[]{
				new LikelihoodColumn(getId())
		};
	}

	private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) {
			super(label);
		}

		public double getDoubleValue() {
			return getLogLikelihood();
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
			Parameter dependentParam = (Parameter) cxo.getChild(Parameter.class);

			cxo = (XMLObject) xo.getChild(INDEPENDENT_VARIABLES);
			Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);

			cxo = (XMLObject) xo.getChild(BASIS_MATRIX);
			DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);

			if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
					(independentParam.getDimension() != designMatrix.getColumnDimension()))
				throw new XMLParseException(
						"dim(" + DEPENDENT_VARIABLES + ") != dim(" + BASIS_MATRIX + " %*% " + INDEPENDENT_VARIABLES + ")"
				);


			String family = xo.getStringAttribute(FAMILY);
			GeneralizedLinearModel glm;
			if (family.compareTo(LOGISTIC_REGRESSION) == 0) {
				glm = new LogisticRegression(dependentParam, independentParam, designMatrix);
			} else
				throw new XMLParseException("Family '" + family + "' is not currently implemented");

			if (glm.requiresScale()) {
				cxo = (XMLObject) xo.getChild(SCALE);
				Parameter scaleParameter = (Parameter) cxo.getChild(Parameter.class);
				if (scaleParameter == null)
					throw new XMLParseException("Family '" + family + "' requires a scale parameter");
				glm.addScaleParameter(scaleParameter);
			}

			return glm;
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
						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				new ElementRule(INDEPENDENT_VARIABLES,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				new ElementRule(BASIS_MATRIX,
						new XMLSyntaxRule[]{new ElementRule(DesignMatrix.class)})
		};

		public String getParserDescription() {
			return "Calculates the generalized linear model likelihood of the dependent parameters given the indepenent parameters and design matrix.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}
	};


}
