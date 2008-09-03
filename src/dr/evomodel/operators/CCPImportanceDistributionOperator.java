/**
 * 
 */
package dr.evomodel.operators;

import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;

import dr.evolution.tree.Clade;
import dr.evolution.tree.ConditionalCladeFrequency;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
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
public class CCPImportanceDistributionOperator extends
		AbstractImportanceDistributionOperator {

	public static final String CCP_IMPORTANCE_DISTRIBUTION_OPERATOR = "CCPImportanceDistributionOperator";

	private int sampleEvery;

	private int samples;

	private int sampleCount;

	private ConditionalCladeFrequency probabilityEstimater;

	private Queue<NodeRef> internalNodes;

	private Map<Integer, NodeRef> externalNodes;

	private boolean burnin = false;

	/**
	 * 
	 */
	public CCPImportanceDistributionOperator(TreeModel tree, double weight,
			int samples, int sampleEvery, double epsilon) {
		super(tree);

		setWeight(weight);
		this.samples = samples;
		this.sampleEvery = sampleEvery;
		probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);

		init();
	}

	/**
	 * 
	 */
	public CCPImportanceDistributionOperator(TreeModel tree, double weight) {
		super(tree);

		setWeight(weight);
		this.samples = 10000;
		sampleEvery = 10;
		probabilityEstimater = new ConditionalCladeFrequency(tree, 1.0);

		init();
	}

	private void init() {
		sampleCount = 0;
		internalNodes = new LinkedList<NodeRef>();
		externalNodes = new HashMap<Integer, NodeRef>();
		fillExternalNodes(tree.getRoot());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.AbstractImportanceSampler#doOperation()
	 */
	@Override
	public double doOperation() throws OperatorFailedException {
		if (!burnin) {
			if (sampleCount < samples * sampleEvery) {
				sampleCount++;
				if (sampleCount % sampleEvery == 0) {
					probabilityEstimater.addTree(tree);
				}
				setAccepted(0);
				setRejected(0);
				setTransitions(0);

				return doUnguidedOperation();

			} else {
				return doImportanceDistributionOperation();
			}
		} else {

			return doUnguidedOperation();

		}
	}

	private double doImportanceDistributionOperation() throws OperatorFailedException {
		final NodeRef root = tree.getRoot();
		BitSet all = new BitSet();
		all.set(0, (tree.getNodeCount() + 1) / 2);
		Clade rootClade = new Clade(all, tree.getNodeHeight(root));

		internalNodes.clear();
		fillInternalNodes(root);
		// remove the root
		internalNodes.poll();

		double prob;
		try {
			prob = createTree(root, rootClade);
		} catch (InvalidTreeException e) {
			throw new OperatorFailedException(e.getMessage());
		}

		// TODO 
		// calculate the probability of the current tree!!!
		// Hr = backward / forward
		
		return prob;
	}

	private void fillInternalNodes(NodeRef node) {
		if (!tree.isExternal(node)) {
			internalNodes.add(node);
			int childCount = tree.getChildCount(node);
			for (int i = 0; i < childCount; i++) {
				fillInternalNodes(tree.getChild(node, i));
			}
		}
	}
	
	private void fillExternalNodes(NodeRef node) {
		if (!tree.isExternal(node)) {
			int childCount = tree.getChildCount(node);
			for (int i = 0; i < childCount; i++) {
				fillExternalNodes(tree.getChild(node, i));
			}
		}
		else {
			String name = tree.getTaxonId(node.getNumber());
			Integer i = probabilityEstimater.getTaxonMap().get(name);
			externalNodes.put(i, node);
		}
	}

	private double createTree(NodeRef node, Clade c) throws InvalidTreeException {
		double prob = 0.0;
		if (c.getSize() == 2) {
			// this clade only contains two tips
			// the split between them is trivial

			int leftTipIndex = c.getBits().nextSetBit(0);
			int rightTipIndex = c.getBits().nextSetBit(leftTipIndex + 1);
			NodeRef leftTip = externalNodes.get(leftTipIndex);
			NodeRef rightTip = externalNodes.get(rightTipIndex);

			tree.beginTreeEdit();
			removeChildren(node);
			tree.addChild(node, leftTip);
			tree.addChild(node, rightTip);
			tree.endTreeEdit();
		} else {
			Clade[] clades = new Clade[2];
			prob = splitClade(c, clades);
			NodeRef leftChild, rightChild;

			if (clades[0].getSize() == 1) {
				int tipIndex = clades[0].getBits().nextSetBit(0);
				leftChild = externalNodes.get(tipIndex);
			} else {
				leftChild = internalNodes.poll();
				// TODO set the node height for the new node
				tree.setNodeHeight(leftChild, tree.getNodeHeight(node) / 2.0);
				prob += createTree(leftChild, clades[0]);
			}

			if (clades[1].getSize() == 1) {
				int tipIndex = clades[1].getBits().nextSetBit(0);
				rightChild = externalNodes.get(tipIndex);
			} else {
				rightChild = internalNodes.poll();
				// TODO set the node height for the new node
				tree.setNodeHeight(rightChild, tree.getNodeHeight(node) / 2.0);
				prob += createTree(rightChild, clades[1]);
			}

			tree.beginTreeEdit();
			removeChildren(node);
			tree.addChild(node, leftChild);
			tree.addChild(node, rightChild);
			tree.endTreeEdit();
		}

		return prob;
	}

	private void removeChildren(NodeRef parent) {
		int childCount = tree.getChildCount(parent);
		for (int i = 0; i < childCount; i++) {
			NodeRef child = tree.getChild(parent, i);
			tree.removeChild(parent, child);
		}
	}

	private double splitClade(Clade c, Clade[] children) {
		return probabilityEstimater.splitClade(c, children);
	}

	/**
	 * @param sampleEvery
	 *            the sampleEvery to set
	 */
	public void setSampleEvery(int sampleEvery) {
		this.sampleEvery = sampleEvery;
	}

	/**
	 * @param samples
	 *            the samples to set
	 */
	public void setSamples(int samples) {
		this.samples = samples;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.AbstractImportanceSampler#getOperatorName()
	 */
	@Override
	public String getOperatorName() {
		return CCP_IMPORTANCE_DISTRIBUTION_OPERATOR;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * dr.inference.operators.AbstractImportanceSampler#getPerformanceSuggestion
	 * ()
	 */
	@Override
	public String getPerformanceSuggestion() {
		if (getAcceptanceProbability() < getMinimumGoodAcceptanceLevel()) {
			return "Try to increase the sample size and/or the steps between each sample.";
		}
		return "";
	}

	public static XMLObjectParser CCP_IMPORTANCE_DISTRIBUTION_OPERATOR_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return CCP_IMPORTANCE_DISTRIBUTION_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			double weight = xo.getDoubleAttribute("weight");
			int samples = xo.getIntegerAttribute("samples");

			double epsilon = 1.0;
			if (xo.hasAttribute("epsilon")) {
				epsilon = xo.getDoubleAttribute("epsilon");
			}

			int sampleEvery = 10;
			if (xo.hasAttribute("sampleEvery")) {
				sampleEvery = xo.getIntegerAttribute("sampleEvery");
			}

			return new CCPImportanceDistributionOperator(treeModel, weight,
					samples, sampleEvery, epsilon);
		}

		//**********************************************************************
		// **
		// AbstractXMLObjectParser implementation
		//**********************************************************************
		// **

		public String getParserDescription() {
			return "This element represents an operator proposing trees from an importance distribution which is created by the conditional clade prbabilities.";
		}

		public Class getReturnType() {
			return CCPImportanceDistributionOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newDoubleRule("weight"),
				AttributeRule.newIntegerRule("samples"),
				AttributeRule.newIntegerRule("sampleEvery", true),
				AttributeRule.newDoubleRule("epsilon", true),
				new ElementRule(TreeModel.class) };

	};
}
