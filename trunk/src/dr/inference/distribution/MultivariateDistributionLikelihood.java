package dr.inference.distribution;

import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.InverseWishartDistribution;
import dr.math.MultivariateDistribution;
import dr.math.MultivariateNormalDistribution;
import dr.math.WishartDistribution;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Marc Suchard
 */
public class MultivariateDistributionLikelihood extends AbstractDistributionLikelihood {

	public static final String MVN_PRIOR = "multivariateNormalPrior";
	public static final String MVN_MEAN = "meanParameter";
	public static final String MVN_PRECISION = "precisionParameter";
	public static final String WISHART_PRIOR = "multivariateWishartPrior";
	public static final String INV_WISHART_PRIOR = "multivariateInverseWishartPrior";
	public static final String DF = "df";
	public static final String SCALE_MATRIX = "scaleMatrix";
	public static final String DATA = "data";

	private MultivariateDistribution distribution;

	public MultivariateDistributionLikelihood(MultivariateDistribution distribution) {
		super(null);
		this.distribution = distribution;
	}

	public void addData(Parameter data) {
		dataList.add(data);
	}

	public String toString() {
		return getId() + "(" + calculateLogLikelihood() + ")";
	}

	protected ArrayList<Parameter> dataList = new ArrayList<Parameter>();


	protected double calculateLogLikelihood() {
		double logL = 0.0;

		for (Parameter parameter : dataList) {
			logL += distribution.logPdf(parameter);
		}
		return logL;
	}

	public MultivariateDistribution getDistribution() {
		return distribution;
	}


	public static XMLObjectParser INV_WISHART_PRIOR_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return INV_WISHART_PRIOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int df = xo.getIntegerAttribute(DF);

			XMLObject cxo = (XMLObject) xo.getChild(SCALE_MATRIX);
			MatrixParameter scaleMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
			InverseWishartDistribution invWishart = new InverseWishartDistribution(df, scaleMatrix.getParameterAsMatrix());

			MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
					invWishart);

			cxo = (XMLObject) xo.getChild(DATA);
			for (int j = 0; j < cxo.getChildCount(); j++) {
				if (cxo.getChild(j) instanceof MatrixParameter) {
					likelihood.addData((MatrixParameter) cxo.getChild(j));
				} else {
					throw new XMLParseException("illegal element in " + xo.getName() + " element " + cxo.getName());
				}
			}


			return likelihood;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(DF),
				new ElementRule(SCALE_MATRIX,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
		};

		public String getParserDescription() {
			return "Calculates the likelihood of some data under an Inverse-Wishart distribution.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}
	};

	public static XMLObjectParser WISHART_PRIOR_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return WISHART_PRIOR;
		}


		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int df = xo.getIntegerAttribute(DF);

			XMLObject cxo = (XMLObject) xo.getChild(SCALE_MATRIX);
			MatrixParameter scaleMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

			MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
					new WishartDistribution(df, scaleMatrix.getParameterAsMatrix())
			);
			cxo = (XMLObject) xo.getChild(DATA);
			for (int j = 0; j < cxo.getChildCount(); j++) {
				if (cxo.getChild(j) instanceof MatrixParameter) {
					likelihood.addData((MatrixParameter) cxo.getChild(j));
					System.err.println("added ");
				} else {
					throw new XMLParseException("illegal element in " + xo.getName() + " element " + cxo.getName());
				}
			}

			return likelihood;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule(DF),
				new ElementRule(SCALE_MATRIX,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
		};

		public String getParserDescription() {
			return "Calculates the likelihood of some data under a Wishart distribution.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}
	};

	public static XMLObjectParser MVN_PRIOR_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MVN_PRIOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//	         double mean = xo.getDoubleAttribute(MVN_MEAN);
//	         double stdev = xo.getDoubleAttribute(STDEV);


			XMLObject cxo = (XMLObject) xo.getChild(MVN_MEAN);
			Parameter mean = (Parameter) cxo.getChild(Parameter.class);

			cxo = (XMLObject) xo.getChild(MVN_PRECISION);
			MatrixParameter precision = (MatrixParameter) cxo.getChild(MatrixParameter.class);

			if (mean.getDimension() != precision.getRowDimension() ||
					mean.getDimension() != precision.getColumnDimension())
				throw new XMLParseException("Mean and precision have wrong dimensions in " + xo.getName() + " element");

			MultivariateDistributionLikelihood likelihood =
					new MultivariateDistributionLikelihood(
							new MultivariateNormalDistribution(mean.getParameterValues(),
									precision.getParameterAsMatrix())
					);
			cxo = (XMLObject) xo.getChild(DATA);
			for (int j = 0; j < cxo.getChildCount(); j++) {
				if (cxo.getChild(j) instanceof Parameter) {
					Parameter data = (Parameter) cxo.getChild(j);
					likelihood.addData(data);
					if (data.getDimension() != mean.getDimension())
						throw new XMLParseException("dim(" + data.getStatisticName() + ") != dim(" + mean.getStatisticName() + ") in " + xo.getName() + "element");
				} else {
					throw new XMLParseException("illegal element in " + xo.getName() + " element");
				}
			}

			return likelihood;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(MVN_MEAN,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				new ElementRule(MVN_PRECISION,
						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
				new ElementRule(DATA,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)})
		};

		public String getParserDescription() {
			return "Calculates the likelihood of some data under a given multivariate-normal distribution.";
		}

		public Class getReturnType() {
			return MultivariateDistributionLikelihood.class;
		}
	};


}
