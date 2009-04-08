package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Hijacks the DiscretizedBranchRates model to allow branch rates to take on any double value
 * This is useful for forming a scaled mixture of normals for the continuous diffusion model
 *
 * @author Marc A. Suchard
 */
public class ArbitraryBranchRates extends DiscretizedBranchRates {


	public static final String ARBITRARY_BRANCH_RATES = "arbitraryBranchRates";
	public static final String RATES = "rates";

	//	private TreeModel tree;
	private Parameter rate;

	public ArbitraryBranchRates(TreeModel tree, Parameter rate) {
		super(tree, rate, null, 1);

//        this.tree = tree;
		this.rate = rate;

		//Force the boundaries of rate
		Parameter.DefaultBounds bound = new Parameter.DefaultBounds(Double.MAX_VALUE, 0, rate.getDimension());
		rate.addBounds(bound);

	}

	protected void setupRates() {
	}

	public double getBranchRate(Tree tree, NodeRef node) {

		int nodeNumber = node.getNumber();
		return rate.getParameterValue(getCategoryIndexFromNodeNumber(nodeNumber));
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return ARBITRARY_BRANCH_RATES;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

			Parameter rateCategoryParameter = (Parameter) xo.getChild(RATES);

			Logger.getLogger("dr.evomodel").info("Using an scaled mixture of normals model.");
//	        Logger.getLogger("dr.evomodel").info("  over sampling = " + overSampling);
//	        Logger.getLogger("dr.evomodel").info("  parametric model = " + distributionModel.getModelName());
			Logger.getLogger("dr.evomodel").info("  rates = " + rateCategoryParameter.getDimension());
			Logger.getLogger("dr.evomodel").info("  NB: Make sure you have a prior on " + rateCategoryParameter.getId() + " and not do use this model in a treeLikelihood");


			return new ArbitraryBranchRates(tree, rateCategoryParameter);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return
					"This element returns an arbitrary rate model." +
							"The branch rates are drawn from an arbitrary distribution determine by the prior.";
		}

		public Class getReturnType() {
			return ArbitraryBranchRates.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//	            AttributeRule.newBooleanRule(SINGLE_ROOT_RATE, true, "Whether only a single rate should be used for the two children branches of the root"),
//	            AttributeRule.newIntegerRule(OVERSAMPLING, true, "The integer factor for oversampling the distribution model (1 means no oversampling)"),
				new ElementRule(TreeModel.class),
//	            new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
				new ElementRule(RATES, Parameter.class, "The rate parameter"),
		};
	};


}
