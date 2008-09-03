/**
 * 
 */
package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.ConditionalCladeFrequency;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
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
 * This class implements a subtree swap operator. The first subtree is chosen
 * randomly and the second one is chosen according to the importance of the new
 * tree. The importance are calculated by the multiplied clade probabilities.
 * 
 */
public class ImportanceSubtreeSwap extends AbstractTreeOperator {

	public static final String IMPORTANCE_SUBTREE_SWAP = "ImportanceSubtreeSwap";

	public final int SAMPLE_EVERY = 10;
	
	private TreeModel tree;

	private int samples;
	
	private int sampleCount = 0;

	private boolean burnin = false;

	private ConditionalCladeFrequency probabilityEstimater;

	/**
	 * 
	 */
	public ImportanceSubtreeSwap(TreeModel tree, double weight, int samples, int epsilon) {
		this.tree = tree;
		setWeight(weight);
		this.samples = samples;
		sampleCount = 0;
		probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
	}
	
	/**
	 * 
	 */
	public ImportanceSubtreeSwap(TreeModel tree, double weight, int samples) {
		this.tree = tree;
		setWeight(weight);
		this.samples = samples;
		sampleCount = 0;
		probabilityEstimater = new ConditionalCladeFrequency(tree, 1.0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
	 */
	@Override
	public double doOperation() throws OperatorFailedException {
		if (!burnin) {
			if (sampleCount < samples * SAMPLE_EVERY) {
				sampleCount++;
				if (sampleCount % SAMPLE_EVERY == 0){
					probabilityEstimater.addTree(tree);					
				}
				setAccepted(0);
				setRejected(0);
				setTransitions(0);
				
				return wideExchange();
				
			} else {
				return importanceExchange();
			}
		} else {
			
			return wideExchange();
			
		 }
	}

	private double wideExchange() throws OperatorFailedException {
		final int nodeCount = tree.getNodeCount();
		final NodeRef root = tree.getRoot();

		NodeRef i;

		do {
			i = tree.getNode(MathUtils.nextInt(nodeCount));
		} while (root == i);

		NodeRef j;
		do {
			j = tree.getNode(MathUtils.nextInt(nodeCount));
		} while (j == i || j == root);

		final NodeRef iP = tree.getParent(i);
		final NodeRef jP = tree.getParent(j);

		if ((iP != jP) && (i != jP) && (j != iP)
				&& (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
				&& (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
			exchangeNodes(tree, i, j, iP, jP);
			return 0.0;
		}

		return 0.0;
	}

	private double importanceExchange()
			throws OperatorFailedException {
		final int nodeCount = tree.getNodeCount();
		final NodeRef root = tree.getRoot();

		NodeRef i;
		int indexI;
		int indexJ;
		
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
		SimpleTree originalTree = new SimpleTree(tree);
		double backward = calculateTreeProbability(originalTree);
		int offset = (int) -backward;
		backward = Math.exp(backward + offset);
		SimpleTree clone;
		for (int n = 0; n < nodeCount; n++) {
			j = tree.getNode(n);
			if (j != root) {
				jP = tree.getParent(j);

				if ((iP != jP) && (i != jP) && (j != iP)
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
						&& (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
					secondNodeIndices.add(n);
					clone = (SimpleTree) originalTree.getCopy();

					swap(clone, clone.getNode(indexI), clone
							.getNode(n));
					double prob = Math.exp(calculateTreeProbability(clone)
							+ offset);
					probabilities.add(prob);
					sum += prob;

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

		// *******************************************
		// assuming we would have chosen j first
		double sumForward2 = 0.0;
		NodeRef k, kP;
		indexJ = secondNodeIndices.get(index);
		for (int n = 0; n < nodeCount; n++) {
			k = tree.getNode(n);
			if (k != root) {
				kP = tree.getParent(k);

				if ((jP != kP) && (j != kP) && (k != jP)
						&& (tree.getNodeHeight(k) < tree.getNodeHeight(jP))
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(kP))) {
					clone = (SimpleTree) originalTree.getCopy();

					swap(clone, clone.getNode(indexJ), clone
							.getNode(n));
					double prob = Math.exp(calculateTreeProbability(clone)
							+ offset);
					sumForward2 += prob;

				}
			}
		}

		exchangeNodes(tree, i, j, iP, jP);
		double forward = probabilities.get(index);

		iP = tree.getParent(i);
		double sumBackward = 0.0;
		for (int n = 0; n < nodeCount; n++) {
			j = tree.getNode(n);
			if (j != root) {
				jP = tree.getParent(j);

				if ((iP != jP) && (i != jP) && (j != iP)
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
						&& (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
					clone = (SimpleTree) originalTree.getCopy();

					swap(clone, clone.getNode(indexI), clone
							.getNode(n));
					double prob = Math.exp(calculateTreeProbability(clone)
							+ offset);
					sumBackward += prob;

				}
			}
		}

		// *******************************************
		// assuming we would have chosen j first
		double sumBackward2 = 0.0;
		j = tree.getNode(secondNodeIndices.get(index));
		jP = tree.getParent(j);
		for (int n = 0; n < nodeCount; n++) {
			k = tree.getNode(n);
			if (k != root) {
				kP = tree.getParent(k);

				if ((jP != kP) && (j != kP) && (k != jP)
						&& (tree.getNodeHeight(k) < tree.getNodeHeight(jP))
						&& (tree.getNodeHeight(j) < tree.getNodeHeight(kP))) {
					clone = (SimpleTree) originalTree.getCopy();

					swap(clone, clone.getNode(indexJ), clone
							.getNode(n));
					double prob = Math.exp(calculateTreeProbability(clone)
							+ offset);
					sumBackward2 += prob;

				}
			}
		}

		double forwardProb = (forward / sum) + (forward / sumForward2);
		double backwardProb = (backward / sumBackward)
				+ (backward / sumBackward2);

		double hastingsRatio = Math.log(backwardProb / forwardProb);

		// throw new OperatorFailedException(
		// "Couldn't find valid wide move on this tree!");

		return hastingsRatio;
	}

	/* exchange subtrees whose root are i and j */
	private SimpleTree swap(SimpleTree tree, NodeRef i, NodeRef j)
			throws OperatorFailedException {

		NodeRef iP = tree.getParent(i);
		NodeRef jP = tree.getParent(j);

		tree.beginTreeEdit();
		tree.removeChild(iP, i);
		tree.removeChild(jP, j);
		tree.addChild(jP, i);
		tree.addChild(iP, j);
		tree.endTreeEdit();

		return tree;
	}

	private double calculateTreeProbability(SimpleTree tree) {
		// return calculateTreeProbabilityMult(tree);
//		return calculateTreeProbabilityLog(tree);
		return probabilityEstimater.getTreeProbability(tree);
//		return 0.0;
	}

	public void setBurnin(boolean burnin) {
		this.burnin = burnin;
	}

	protected double evaluate(Likelihood likelihood, Prior prior) {

		double logPosterior = 0.0;

		if (prior != null) {
			final double logPrior = prior.getLogPrior(likelihood.getModel());

			if (logPrior == Double.NEGATIVE_INFINITY) {
				return Double.NEGATIVE_INFINITY;
			}

			logPosterior += logPrior;
		}

		final double logLikelihood = likelihood.getLogLikelihood();

		if (Double.isNaN(logLikelihood)) {
			return Double.NEGATIVE_INFINITY;
		}
		// System.err.println("** " + logPosterior + " + " + logLikelihood + " =
		// " + (logPosterior + logLikelihood));
		logPosterior += logLikelihood;

		return logPosterior;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
	 */
	@Override
	public String getOperatorName() {
		return IMPORTANCE_SUBTREE_SWAP;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
	 */
	public String getPerformanceSuggestion() {
		// TODO Auto-generated method stub
		return "";
	}

	public static XMLObjectParser IMPORTANCE_SUBTREE_SWAP_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return IMPORTANCE_SUBTREE_SWAP;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			double weight = xo.getDoubleAttribute("weight");
			int samples = xo.getIntegerAttribute("samples");

			return new ImportanceSubtreeSwap(treeModel, weight, samples);
		}

		// ************************************************************************
		// AbstractXMLObjectParser implementation
		// ************************************************************************

		public String getParserDescription() {
			return "This element represents a importance guided subtree swap operator. "
					+ "This operator swaps a random subtree with a second subtree guided by an importance distribution.";
		}

		public Class getReturnType() {
			return ImportanceSubtreeSwap.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newDoubleRule("weight"),
				AttributeRule.newIntegerRule("samples"),
				new ElementRule(TreeModel.class) };

	};

}
