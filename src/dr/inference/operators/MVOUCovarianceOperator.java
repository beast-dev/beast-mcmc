package dr.inference.operators;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.WishartDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class MVOUCovarianceOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

	public static final String MVOU_OPERATOR = "mvouOperator";
	public static final String MIXING_FACTOR = "mixingFactor";
	public static final String VARIANCE_MATRIX = "varMatrix";
	public static final String PRIOR_DF = "priorDf";

	private int mode;
	private double mixingFactor;
	private MatrixParameter varMatrix;
	private int dim;

	private MatrixParameter precisionParam;
	private WishartDistribution priorDistribution;
	private int priorDf;
	private double[][] I;
	private Matrix Iinv;


	public MVOUCovarianceOperator(double mixingFactor,
	                              MatrixParameter varMatrix,
	                              int priorDf,
	                              double weight, int mode) {
		super();
		this.mode = mode;
		this.mixingFactor = mixingFactor;
		this.varMatrix = varMatrix;
		this.priorDf = priorDf;
		setWeight(weight);
		dim = varMatrix.getColumnDimension();
		I = new double[dim][dim];
		for (int i = 0; i < dim; i++)
			I[i][i] = 1.0;
		Iinv = new Matrix(I).inverse();
	}

	public double doOperation() throws OperatorFailedException {

		double[][] draw = WishartDistribution.nextWishart(priorDf, I);
		double[][] oldValue = varMatrix.getParameterAsMatrix();
		for (int i = 0; i < dim; i++) {
			Parameter column = varMatrix.getParameter(i);
			for (int j = 0; j < dim; j++)
				column.setParameterValue(j,
						mixingFactor *
								oldValue[j][i]
								+ (1.0 - mixingFactor) * draw[j][i]
				);

		}
//        varMatrix.fireParameterChangedEvent();
		// calculate Hastings ratio
		Matrix forwardDrawMatrix = new Matrix(draw);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				oldValue[i][j] -= mixingFactor * varMatrix.getParameterValue(i, j);
				oldValue[i][j] /= 1.0 - mixingFactor;

			}
		}
		Matrix backwardDrawMatrix = new Matrix(oldValue);

		return WishartDistribution.logPdf(backwardDrawMatrix, Iinv, priorDf, dim, 0)
				- WishartDistribution.logPdf(forwardDrawMatrix, Iinv, priorDf, dim, 0);
	}

	//MCMCOperator INTERFACE
	public final String getOperatorName() {
		return MVOU_OPERATOR + "(" +
				varMatrix.getId() + ")";
	}

	public double getCoercableParameter() {
		return Math.log(mixingFactor / (1.0 - mixingFactor));
	}

	public void setCoercableParameter(double value) {
		mixingFactor = Math.exp(value) / (1.0 + Math.exp(value));
	}

	public double getRawParameter() {
		return mixingFactor;
	}

	public int getMode() {
		return mode;
	}

	public double getMixingFactor() {
		return mixingFactor;
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

	public final String getPerformanceSuggestion() {

		double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
		double targetProb = getTargetAcceptanceProbability();
		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
		double sf = OperatorUtils.optimizeWindowSize(mixingFactor, prob, targetProb);
		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try setting mixingFactor to about " + formatter.format(sf);
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try setting mixingFactor to about " + formatter.format(sf);
		} else return "";
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return MVOU_OPERATOR;
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

			double weight = xo.getDoubleAttribute(WEIGHT);
			double mixingFactor = xo.getDoubleAttribute(MIXING_FACTOR);
			int priorDf = xo.getIntegerAttribute(PRIOR_DF);

			if (mixingFactor <= 0.0 || mixingFactor >= 1.0) {
				throw new XMLParseException("mixingFactor must be greater than 0.0 and less thatn 1.0");
			}

//            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

//            XMLObject cxo = (XMLObject) xo.getChild(VARIANCE_MATRIX);
			MatrixParameter varMatrix = (MatrixParameter) xo.getChild(MatrixParameter.class);

			// Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

			if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
				throw new XMLParseException("The variance matrix is not square");

//            if (varMatrix.getColumnDimension() != parameter.getDimension())
//                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

			return new MVOUCovarianceOperator(mixingFactor, varMatrix, priorDf, weight, mode);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns junk.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(MIXING_FACTOR),
				AttributeRule.newIntegerRule(PRIOR_DF),
				AttributeRule.newDoubleRule(WEIGHT),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//                new ElementRule(Parameter.class),
//                new ElementRule(VARIANCE_MATRIX,
//                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),

				new ElementRule(MatrixParameter.class)

		};

	};
}
