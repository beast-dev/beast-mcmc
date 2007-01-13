package dr.inference.operators;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.xml.*;
//import numericalMethods.algebra.linear.decompose.Cholesky;


/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Dec 29, 2006
 * Time: 3:49:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultivariateNormalOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

	public static final String MVN_OPERATOR = "mvnOperator";
	public static final String SCALE_FACTOR = "scaleFactor";
	public static final String VARIANCE_MATRIX = "varMatrix";

	private int mode;
	private double scaleFactor;
	private Parameter parameter;
	private int dim;

	private double[][] cholesky;

	public MultivariateNormalOperator(Parameter parameter, double scaleFactor,
	                                  MatrixParameter varMatrix, int weight, int mode) {
		super();
		this.mode = mode;
		this.scaleFactor = scaleFactor;
		this.parameter = parameter;
		this.weight = weight;
		dim = parameter.getDimension();
		cholesky = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++)
				cholesky[i][j] = cholesky[j][i] = varMatrix.getParameterValue(i, j);
		}
		//Cholesky.decompose(cholesky, cholesky);

		try {
			cholesky = (new CholeskyDecomposition(cholesky)).getL();
		} catch (IllegalDimension illegalDimension) {
			// todo - check for square variance matrix before here
		}

	}

	public double doOperation() throws OperatorFailedException {

		double[] x = parameter.getParameterValues();
		double[] epsilon = new double[dim];
		//double[] y = new double[dim];
		for (int i = 0; i < dim; i++)
			epsilon[i] = scaleFactor * MathUtils.nextGaussian();

		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++) {
				x[i] += cholesky[i][j] * epsilon[j];
			}
			parameter.setParameterValue(i, x[i]);
//            System.out.println(i+" : "+x[i]);
		}
//                    System.exit(-1);
		return 0;
	}

	//MCMCOperator INTERFACE
	public final String getOperatorName() {
		return parameter.getParameterName();
	}

	public double getCoercableParameter() {
		return Math.log(scaleFactor);
	}

	public void setCoercableParameter(double value) {
		scaleFactor = Math.exp(value);
	}

	public double getRawParameter() {
		return scaleFactor;
	}

	public int getMode() {
		return mode;
	}

	public double getScaleFactor() {
		return scaleFactor;
	}

	public double getTargetAcceptanceProbability() {
		return 0.234;
	}

	public double getMinimumAcceptanceLevel() {
		return 0.1;
	}

	public double getMaximumAcceptanceLevel() {
		return 0.4;
	}

	public double getMinimumGoodAcceptanceLevel() {
		return 0.20;
	}

	public double getMaximumGoodAcceptanceLevel() {
		return 0.30;
	}

	public int getWeight() {
		return weight;
	}

	public void setWeight(int w) {
		weight = w;
	}

	public final String getPerformanceSuggestion() {

		double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
		double targetProb = getTargetAcceptanceProbability();
		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
		double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);
		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else return "";
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return MVN_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int mode = CoercableMCMCOperator.DEFAULT;

			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
					mode = CoercableMCMCOperator.COERCION_ON;
				} else {
					mode = CoercableMCMCOperator.COERCION_OFF;
				}
			}

			int weight = xo.getIntegerAttribute(WEIGHT);
			double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

			if (scaleFactor <= 0.0) {
				throw new XMLParseException("scaleFactor must be greater than 0.0");
			}

			Parameter parameter = (Parameter) xo.getChild(Parameter.class);


			XMLObject cxo = (XMLObject) xo.getChild(VARIANCE_MATRIX);
			MatrixParameter varMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

			// Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

			if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
				throw new XMLParseException("The variance matrix is not square");

			if (varMatrix.getColumnDimension() != parameter.getDimension())
				throw new XMLParseException("The parameter and variance matrix have differing dimensions");

			return new MultivariateNormalOperator(parameter, scaleFactor, varMatrix, weight, mode);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a multivariate normal random walk operator on a given parameter.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(SCALE_FACTOR),
				AttributeRule.newIntegerRule(WEIGHT),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
				new ElementRule(Parameter.class),
				new ElementRule(VARIANCE_MATRIX,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)})

		};

	};
}
