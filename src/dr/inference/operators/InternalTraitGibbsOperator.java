package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.math.MathUtils;
import dr.math.MultivariateNormalDistribution;
import dr.xml.*;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jun 13, 2007
 * Time: 9:18:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class InternalTraitGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

	public static final String GIBBS_OPERATOR = "internalTraitGibbsOperator";
//	public static final String PRECISION_MATRIX = "precisionMatrix";
//	public static final String TREE_MODEL = "treeModel";
//	public static final String OUTCOME  = "outcome";
//	public static final String MEAN = "mean";
//	public static final String PRIOR = "prior";

//	private Parameter outcomeParam;
//	private Parameter meanParam;
//	private MatrixParameter precisionParam;
//	private WishartDistribution priorDistribution;
	//	private int priorDf;
	//	private double[][] priorScaleMatrix;
	//	private SymmetricMatrix priorScaleMatrix;
	private TreeModel treeModel;
	private MultivariateDiffusionModel diffusionModel;
	private MatrixParameter precisionMatrixParameter;
	private int dim;
//	private int numberObservations;
//	private int weight;

	public InternalTraitGibbsOperator(TreeModel treeModel, MultivariateDiffusionModel diffusionModel) {
		super();
		this.treeModel = treeModel;
		this.diffusionModel = diffusionModel;
//		System.err.println("D: "+diffusionModel.toString());
//		System.exit(0);
		precisionMatrixParameter = diffusionModel.getPrecisionMatrixParameter();
		dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), "trait").length;
		// todo need to fix trait name
	}


	public int getStepCount() {
		return 1;
	}


	public double doOperation() throws OperatorFailedException {

		double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

		NodeRef node = treeModel.getRoot();
//		System.err.println("CNT: "+treeModel.getInternalNodeCount());
		while (node == treeModel.getRoot()) {
			node = treeModel.getInternalNode(MathUtils.nextInt(
					treeModel.getInternalNodeCount()));
		} // select any internal node but the root.

		NodeRef parent = treeModel.getParent(node);
//		int childCount = treeModel.getChildCount(node);

		double[] weightedAverage = new double[dim];
//		double weightTotal = 0;

//		double[] weights = new double[childCount+1];
		double weight = 1.0 / treeModel.getBranchLength(node);
		double[] trait = treeModel.getMultivariateNodeTrait(parent, "trait");

		for (int i = 0; i < dim; i++)
			weightedAverage[i] = trait[i] * weight;
//		weightedAverage[i] = i;

		double weightTotal = weight;
		for (int j = 0; j < treeModel.getChildCount(node); j++) {
			NodeRef child = treeModel.getChild(node, j);
			trait = treeModel.getMultivariateNodeTrait(child, "trait");
			weight = 1.0 / treeModel.getBranchLength(child);
			for (int i = 0; i < dim; i++)
				weightedAverage[i] += trait[i] * weight;
			weightTotal += weight;
		}

		for (int i = 0; i < dim; i++) {
			weightedAverage[i] /= weightTotal;
			for (int j = i; j < dim; j++)
				precision[j][i] = precision[i][j] * weightTotal;
		}

		double[] draw = MultivariateNormalDistribution.nextMultivariateNormal(
				weightedAverage, precision);
//	    treeModel.setMultivariateTrait(node,weightedAverage);
		treeModel.setMultivariateTrait(node, "trait", draw);
//		for(int i=0; i<dim; i++)

		return 0;  //To change body of implemented methods use File | Settings | File Templates.

	}

	public String getPerformanceSuggestion() {
		return null;  //To change body of implemented methods use File | Settings | File Templates.
	}

	public String getOperatorName() {
		return GIBBS_OPERATOR;
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return GIBBS_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int weight = xo.getIntegerAttribute(WEIGHT);

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

			MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);

			return new InternalTraitGibbsOperator(treeModel, diffusionModel);
//					precMatrix,(WishartDistribution)prior.getDistribution(),treeModel, weight);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a multivariate Gibbs operator on an internal node trait.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//				AttributeRule.newDoubleRule(SCALE_FACTOR),
				AttributeRule.newIntegerRule(WEIGHT),
				new ElementRule(TreeModel.class),
				new ElementRule(MultivariateDiffusionModel.class)
		};

	};

}
