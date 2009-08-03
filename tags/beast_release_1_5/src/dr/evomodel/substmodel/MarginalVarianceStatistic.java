package dr.evomodel.substmodel;

import dr.evomodel.substmodel.MultivariateOUModel;
import dr.xml.*;
import dr.inference.model.Statistic;

/**
 * @author Marc A. Suchard
 */
public class MarginalVarianceStatistic extends Statistic.Abstract {

	public static final String VARIANCE_STATISTIC = "marginalVariance";

	public MarginalVarianceStatistic(MultivariateOUModel mvou) {
		this.mvou = mvou;
	}

	public int getDimension() {
		return mvou.getDimension();
	}

	public double getStatisticValue(int dim) {
		return mvou.getStatisticValue(dim);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return VARIANCE_STATISTIC;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			MultivariateOUModel mvou = (MultivariateOUModel) xo.getChild(MultivariateOUModel.class);
			return new MarginalVarianceStatistic(mvou);

		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the matrix inverse of the child statistic.";
		}

		public Class getReturnType() {
			return MarginalVarianceStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(MultivariateOUModel.class)
		};
	};

	private MultivariateOUModel mvou;
}
