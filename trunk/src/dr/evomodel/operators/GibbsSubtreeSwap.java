/**
 * 
 */
package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.operators.OperatorFailedException;
import dr.inference.prior.Prior;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Sebastian Hoehna
 * 
 */
public class GibbsSubtreeSwap extends SimpleGibbsOperator {

	public static final String GIBBS_SUBTREE_EXCHANGE = "GibbsSubtreeExchange";

	private TreeModel tree;
	
	private Prior prior;
	
	private Likelihood likelihood;

	/**
	 * 
	 */
	public GibbsSubtreeSwap(TreeModel tree, double weight) {
		this.tree = tree;
		setWeight(weight);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.evomodel.operators.SimpleGibbsOperator#getStepCount()
	 */
	@Override
	public int getStepCount() {
		return 1;
	}

	public double doOperation(Prior prior,
            Likelihood likelihood) throws OperatorFailedException {

		this.prior = prior;
		this.likelihood = likelihood;
		
		int tipCount = tree.getExternalNodeCount();

		try {
			wide();
		} catch (InvalidTreeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		if (tree.getExternalNodeCount() != tipCount) {
			throw new RuntimeException("Lost some tips in "
					+ "WIDE mode.");
		}

		return 0.0;
	}

//	private double forward;
//	private double sum = 0.0;
//	private double backward;
//	private double sumBackward = 0.0;
//	private double forwardProb;
//	private double backwardProb;
//	private double hastingsRatio;
//	private double forwardLikelihood;
//	private double backwardLikelihood;

	/**
	 * WARNING: Assumes strictly bifurcating tree.
	 * @throws InvalidTreeException 
	 */
	public void wide() throws OperatorFailedException, InvalidTreeException {

		final int nodeCount = tree.getNodeCount();
		final NodeRef root = tree.getRoot();

		NodeRef i;
		int indexI;

		do {
			indexI = MathUtils.nextInt(nodeCount);
			i = tree.getNode(indexI);
		} while (root == i
				|| (tree.getParent(i) == root && tree.getNodeHeight(i) > tree
						.getNodeHeight(getOtherChild(tree, tree.getParent(i), i))));

		List<Integer> secondNodeIndices = new ArrayList<Integer>();
		List<Double> probabilities = new ArrayList<Double>();
		NodeRef j, iP, jP;
		iP = tree.getParent(i);
		double sum = 0.0;
		double backwardLikelihood = calculateTreeLikelihood(prior, likelihood, tree);
		int offset = (int) -backwardLikelihood;
		double backward = Math.exp(backwardLikelihood + offset);
		for (int n = 0; n < nodeCount; n++) {
			j = tree.getNode(n);
			if (j != root) {
				jP = tree.getParent(j);

//				if ((iP != jP) && (i != jP) && (j != iP)
				if ((i != j) && (i != jP) && (j != iP)
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
						&& (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
					secondNodeIndices.add(n);

					swap(tree, tree.getNode(indexI), tree.getNode(n));
					double prob = Math.exp(calculateTreeLikelihood(prior,
							likelihood, tree)
							+ offset);
					probabilities.add(prob);
					sum += prob;

					undoSwap(tree, tree.getNode(indexI), tree.getNode(n));
				}
			}
		}

		double ran = Math.random() * sum;
		int index = 0;
		while (ran > 0.0) {
			ran -= probabilities.get(index);
			index++;
		}
		index--;

		j = tree.getNode(secondNodeIndices.get(index));
		jP = tree.getParent(j);

		exchangeNodes(tree, i, j, iP, jP);
		double forward = probabilities.get(index);

		iP = tree.getParent(i);
		double sumBackward = 0.0;
		for (int n = 0; n < nodeCount; n++) {
			j = tree.getNode(n);
			if (j != root) {
				jP = tree.getParent(j);

//				if ((iP != jP) && (i != jP) && (j != iP)
				if ((i != j) && (i != jP) && (j != iP)
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
						&& (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {

					swap(tree, tree.getNode(indexI), tree.getNode(n));
					double prob = Math.exp(calculateTreeLikelihood(prior,
							likelihood, tree)
							+ offset);
					sumBackward += prob;

					undoSwap(tree, tree.getNode(indexI), tree.getNode(n));
				}
			}
		}
		
		tree.pushTreeChangedEvent(iP);
		tree.pushTreeChangedEvent(jP);
		

		double forwardProb = (forward / sum);
		double backwardProb = (backward / sumBackward);
		
		double hastingsRatio = Math.log(backwardProb / forwardProb);

		//return hastingsRatio;
	}	
	
	private double calculateTreeLikelihood(Prior prior, Likelihood likelihood,
			TreeModel tree) {
		return evaluate(likelihood, prior);
	}
	
	/**
	 * @param tree   the tree
	 * @param parent the parent
	 * @param child  the child that you want the sister of
	 * @return the other child of the given parent.
	 */
	protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

	    if (tree.getChild(parent, 0) == child) {
	        return tree.getChild(parent, 1);
	    } else {
	        return tree.getChild(parent, 0);
	    }
	}
	
	/* exchange subtrees whose root are i and j */
	private TreeModel swap(TreeModel tree, NodeRef i, NodeRef j)
			throws OperatorFailedException, InvalidTreeException {

		NodeRef iP = tree.getParent(i);
		NodeRef jP = tree.getParent(j);

		tree.beginTreeEdit();
		tree.removeChild(iP, i);
		tree.removeChild(jP, j);
		tree.addChild(jP, i);
		tree.addChild(iP, j);

		tree.pushTreeChangedEvent(iP);
		tree.pushTreeChangedEvent(jP);
		
		tree.endTreeEdit();

		return tree;
	}

	/* undo the subtree exchange */
	private TreeModel undoSwap(TreeModel tree, NodeRef i, NodeRef j)
			throws OperatorFailedException, InvalidTreeException {
		return swap(tree, i, j);
	}
	
	/* exchange subtrees whose root are i and j */
	protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
	                             NodeRef iP, NodeRef jP) throws OperatorFailedException {

	    tree.beginTreeEdit();
	    tree.removeChild(iP, i);
	    tree.removeChild(jP, j);
	    tree.addChild(jP, i);
	    tree.addChild(iP, j);

	    try {
	        tree.endTreeEdit();
	    } catch (MutableTree.InvalidTreeException ite) {
	        throw new OperatorFailedException(ite.toString());
	    }
	}

	public String getOperatorName() {
		return "Gibbs Subtree Exchange";
	}

	public static XMLObjectParser GIBBS_SUBTREE_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return GIBBS_SUBTREE_EXCHANGE;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			double weight = xo.getDoubleAttribute("weight");

			return new GibbsSubtreeSwap(treeModel, weight);
		}

		// ************************************************************************
		// AbstractXMLObjectParser implementation
		// ************************************************************************

		public String getParserDescription() {
			return "This element represents a Gibbs wide exchange operator. "
					+ "This operator swaps two subtrees chosen to their posterior probaility.";
		}

		public Class getReturnType() {
			return GibbsSubtreeSwap.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newDoubleRule("weight"),
				new ElementRule(TreeModel.class) };

	};

}
