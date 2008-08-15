/**
 * 
 */
package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.ConditionalCladeFrequency;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author shhn001
 *
 */
public class ImportancePruneAndRegraft extends AbstractTreeOperator {

	public static final String IMPORTANCE_PRUNE_AND_REGRAFT = "ImportancePruneAndRegraft";

	private TreeModel tree;

	private int samples;

	private int sampleCount;

	private boolean burnin = false;
	
	private ConditionalCladeFrequency probabilityEstimater;

	/**
	 * 
	 */
	public ImportancePruneAndRegraft(TreeModel tree, double weight, int samples, int epsilon) {
		this.tree = tree;
		setWeight(weight);
		this.samples = samples;
		sampleCount = 0;
		probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
	}
	
	/**
	 * 
	 */
	public ImportancePruneAndRegraft(TreeModel tree, double weight, int samples) {
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
			if (sampleCount < samples) {
				probabilityEstimater.addTree(tree);
				setAccepted(0);
				setRejected(0);
				
				return fixedNodePruneAndRegraft();
				
			} else {
				return importancePruneAndRegraft();
			}
		} else {
			
			return fixedNodePruneAndRegraft();
			
		 }
	}

	private double fixedNodePruneAndRegraft() throws OperatorFailedException {
		NodeRef iGrandfather, iBrother;
	      double heightFather;
	      final int tipCount = tree.getExternalNodeCount();

	      final int nNodes = tree.getNodeCount();
	      final NodeRef root = tree.getRoot();

	      NodeRef i;

	      int MAX_TRIES = 1000;

	      for (int tries = 0; tries < MAX_TRIES; ++tries) {
	         // get a random node whose father is not the root - otherwise
	         // the operation is not possible
	         do {
	            i = tree.getNode(MathUtils.nextInt(nNodes));
	         } while (root == i || tree.getParent(i) == root);

	         // int childIndex = (MathUtils.nextDouble() >= 0.5 ? 1 : 0);
	         // int otherChildIndex = 1 - childIndex;
	         // NodeRef iOtherChild = tree.getChild(i, otherChildIndex);

	         NodeRef iFather = tree.getParent(i);
	         iGrandfather = tree.getParent(iFather);
	         iBrother = getOtherChild(tree, iFather, i);
	         heightFather = tree.getNodeHeight(iFather);

	         // NodeRef newChild = getRandomNode(possibleChilds, iFather);
	         NodeRef newChild = tree.getNode(MathUtils.nextInt(nNodes));

	         if (tree.getNodeHeight(newChild) < heightFather
	               && root != newChild
	               && tree.getNodeHeight(tree.getParent(newChild)) > heightFather
	               && newChild != iFather
	               && tree.getParent(newChild) != iFather) {
	            NodeRef newGrandfather = tree.getParent(newChild);

	            tree.beginTreeEdit();

	            // prune
	            tree.removeChild(iFather, iBrother);
	            tree.removeChild(iGrandfather, iFather);
	            tree.addChild(iGrandfather, iBrother);

	            // reattach
	            tree.removeChild(newGrandfather, newChild);
	            tree.addChild(iFather, newChild);
	            tree.addChild(newGrandfather, iFather);

	            // ****************************************************

	            try {
	               tree.endTreeEdit();
	            } catch (MutableTree.InvalidTreeException ite) {
	               throw new OperatorFailedException(ite.toString());
	            }

	            tree.pushTreeChangedEvent(i);
	            tree.pushTreeChangedEvent(iBrother);

	            assert tree.getExternalNodeCount() == tipCount;
	            
	            return 0.0;
	         }
	      }

	      throw new OperatorFailedException("Couldn't find valid SPR move on this tree!");
	}

	public double oldScore;
	public double oldScore2;
	public double adjustedLogScore;

	private double importancePruneAndRegraft()
			throws OperatorFailedException {
		tree.storeModelState();
		
		final int nodeCount = tree.getNodeCount();
		final NodeRef root = tree.getRoot();		

		NodeRef i;
		int indexI;

		do {
			indexI = MathUtils.nextInt(nodeCount);
			i = tree.getNode(indexI);
		} while (root == i || tree.getParent(i) == root);

		List<Integer> secondNodeIndices = new ArrayList<Integer>();
		List<Double> probabilities = new ArrayList<Double>();
		NodeRef j, iP, jP;
		iP = tree.getParent(i);
		double sum = 0.0;
				
		SimpleTree originalTree = new SimpleTree(tree);
		double backwardLikelihood = calculateTreeProbability(originalTree);
		int offset = (int) -backwardLikelihood;
		double backward = Math.exp(backwardLikelihood + offset);
		NodeRef oldBrother;
		oldBrother = getOtherChild(tree, iP, i);
		SimpleTree clone;
		for (int n = 0; n < nodeCount; n++) {
			j = tree.getNode(n);
			if (j != root) {
				jP = tree.getParent(j);

				if (j == oldBrother){					
					secondNodeIndices.add(n);
					probabilities.add(backward);
					sum += backward;
				} else if ((i != j) && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
						&& (tree.getNodeHeight(iP) < tree.getNodeHeight(jP))) {
					secondNodeIndices.add(n);
					clone = (SimpleTree) originalTree.getCopy();

					SimpleNode sI = (SimpleNode)clone.getNode(indexI);
					SimpleNode sIP = (SimpleNode)clone.getParent(sI);
					SimpleNode sJ = (SimpleNode)clone.getNode(n);
					SimpleNode sJP = (SimpleNode)clone.getParent(sJ);
					pruneAndRegraft(clone, sI, sIP, sJ, sJP);
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

		if (iP != jP){
			pruneAndRegraft(tree, i, iP, j, jP);
		}
		double forward = probabilities.get(index);		

		double forwardProb = (forward / sum);
		double backwardProb = (backward / sum);
		double hastingsRatio = Math.log(backwardProb / forwardProb);

		return hastingsRatio;
	}
	
	private void pruneAndRegraft(TreeModel tree, NodeRef i, NodeRef iP, NodeRef j, NodeRef jP) throws OperatorFailedException{
		tree.beginTreeEdit();

		// the grandfather
		NodeRef iG = tree.getParent(iP);
		// the brother
		NodeRef iB = getOtherChild(tree, iP, i);
        // prune
        tree.removeChild(iP, iB);
        tree.removeChild(iG, iP);
        tree.addChild(iG, iB);

        // reattach
        tree.removeChild(jP, j);
        tree.addChild(iP, j);
        tree.addChild(jP, iP);

        // ****************************************************

        try {
           tree.endTreeEdit();
        } catch (MutableTree.InvalidTreeException ite) {
           throw new OperatorFailedException(ite.toString());
        }

        tree.pushTreeChangedEvent(i);
	}
	
	private void pruneAndRegraft(SimpleTree tree, SimpleNode i, SimpleNode iP, SimpleNode j, SimpleNode jP) throws OperatorFailedException{
		tree.beginTreeEdit();

		// the grandfather
		NodeRef iG = tree.getParent(iP);
		// the brother
		NodeRef iB = getOtherChild(tree, iP, i);
        // prune
        tree.removeChild(iP, iB);
        tree.removeChild(iG, iP);
        tree.addChild(iG, iB);

        // reattach
        tree.removeChild(jP, j);
        tree.addChild(iP, j);
        tree.addChild(jP, iP);

        // ****************************************************

        tree.endTreeEdit();
	}
	
	private double calculateTreeProbability(SimpleTree tree) {
		// return calculateTreeProbabilityMult(tree);
//		return calculateTreeProbabilityLog(tree);
		return probabilityEstimater.getTreeProbability(tree);
//		return 10.5;
	}

	public void setBurnin(boolean burnin) {
		this.burnin = burnin;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
	 */
	@Override
	public String getOperatorName() {
		return IMPORTANCE_PRUNE_AND_REGRAFT;
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

	public static XMLObjectParser IMPORTANCE_PRUNE_AND_REGRAFT_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return IMPORTANCE_PRUNE_AND_REGRAFT;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			double weight = xo.getDoubleAttribute("weight");
			int samples = xo.getIntegerAttribute("samples");

			return new ImportancePruneAndRegraft(treeModel, weight, samples);
		}

		// ************************************************************************
		// AbstractXMLObjectParser implementation
		// ************************************************************************

		public String getParserDescription() {
			return "This element represents a importance guided prune and regraft operator. "
					+ "This operator prunes a random subtree and regrafts it below a node chosen by an importance distribution.";
		}

		public Class getReturnType() {
			return ImportancePruneAndRegraft.class;
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
