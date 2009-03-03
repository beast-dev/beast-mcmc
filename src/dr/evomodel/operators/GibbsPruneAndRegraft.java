/**
 * 
 */
package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
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
public class GibbsPruneAndRegraft extends SimpleGibbsOperator {

    public static final String GIBBS_PRUNE_AND_REGRAFT = "GibbsPruneAndRegraft";

    private int MAX_DISTANCE = 10;

    private TreeModel tree;

    private int[] distances;

    private double[] scores;

    private boolean pruned = true;

    /**
	 * 
	 */
    public GibbsPruneAndRegraft(TreeModel tree, boolean pruned, double weight) {
	this.tree = tree;
	this.pruned = pruned;
	setWeight(weight);
	scores = new double[tree.getNodeCount()];
	MAX_DISTANCE = tree.getNodeCount() / 10;
    }

    /*
     * (non-Javadoc)
     * 
     * @see
     * dr.evomodel.operators.SimpleGibbsOperator#doOperation(dr.inference.prior
     * .Prior, dr.inference.model.Likelihood)
     */
    @Override
    public double doOperation(Prior prior, Likelihood likelihood)
	    throws OperatorFailedException {
	if (pruned) {
	    return prunedGibbsProposal(prior, likelihood);
	} else {
	    return GibbsProposal(prior, likelihood);
	}
    }

    private double GibbsProposal(Prior prior, Likelihood likelihood)
	    throws OperatorFailedException {
	final int nodeCount = tree.getNodeCount();
	final NodeRef root = tree.getRoot();

	NodeRef i;

	do {
	    int indexI = MathUtils.nextInt(nodeCount);
	    i = tree.getNode(indexI);
	} while (root == i || tree.getParent(i) == root);

	List<Integer> secondNodeIndices = new ArrayList<Integer>();
	List<Double> probabilities = new ArrayList<Double>();
	NodeRef j, jP;
	final NodeRef iP = tree.getParent(i);
	final double heightIP = tree.getNodeHeight(iP);
	double sum = 0.0;
	double backwardLikelihood = calculateTreeLikelihood(prior, likelihood,
		tree);
	int offset = (int) -backwardLikelihood;
	double backward = Math.exp(backwardLikelihood + offset);
	final NodeRef oldBrother = getOtherChild(tree, iP, i);
	final NodeRef oldGrandfather = tree.getParent(iP);
	for (int n = 0; n < nodeCount; n++) {
	    j = tree.getNode(n);
	    if (j != root) {
		jP = tree.getParent(j);

		if (j == oldBrother) {
		    // secondNodeIndices.add(n);
		    // probabilities.add(backward);
		    // sum += backward;
		} else if ((i != j) && (tree.getNodeHeight(j) < heightIP)
			&& (heightIP < tree.getNodeHeight(jP))) {
		    secondNodeIndices.add(n);

		    pruneAndRegraft(tree, i, iP, j, jP);
		    double prob = Math.exp(calculateTreeLikelihood(prior,
			    likelihood, tree)
			    + offset);
		    probabilities.add(prob);
		    sum += prob;

		    pruneAndRegraft(tree, i, iP, oldBrother, oldGrandfather);
		}
	    }
	}

	if (sum <= 1E-100) {
	    // hack
	    // the proposals have such a small likelihood that they can be
	    // neglected
	    throw new OperatorFailedException(
		    "Couldn't find another proposal with a decent likelihood.");
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

	// int distance = getNodeDistance(iP, jP);
	// distances[distance]++;

	pruneAndRegraft(tree, i, iP, j, jP);

	double forward = probabilities.get(index);

	// tree.pushTreeChangedEvent(jP);
	// tree.pushTreeChangedEvent(oldGrandfather);
	tree.pushTreeChangedEvent(i);

	double forwardProb = (forward / sum);
	double backwardProb = (backward / (sum - forward + backward));
	double hastingsRatio = Math.log(backwardProb / forwardProb);

	return hastingsRatio;
    }

    private double prunedGibbsProposal(Prior prior, Likelihood likelihood)
	    throws OperatorFailedException {
	final int nodeCount = tree.getNodeCount();
	final NodeRef root = tree.getRoot();

	for (int i = 0; i < nodeCount; i++) {
	    scores[i] = Double.NEGATIVE_INFINITY;
	}

	NodeRef i;

	do {
	    int indexI = MathUtils.nextInt(nodeCount);
	    i = tree.getNode(indexI);
	} while (root == i || tree.getParent(i) == root);

	List<Integer> secondNodeIndices = new ArrayList<Integer>();
	List<Double> probabilities = new ArrayList<Double>();
	NodeRef j, jP;
	final NodeRef iP = tree.getParent(i);
	final double heightIP = tree.getNodeHeight(iP);
	double sum = 0.0;
	double backwardLikelihood = calculateTreeLikelihood(prior, likelihood,
		tree);
	int offset = (int) -backwardLikelihood;
	double backward = Math.exp(backwardLikelihood + offset);
	final NodeRef oldBrother = getOtherChild(tree, iP, i);
	final NodeRef oldGrandfather = tree.getParent(iP);
	for (int n = 0; n < nodeCount; n++) {
	    j = tree.getNode(n);
	    if (j != root) {
		jP = tree.getParent(j);

		if (j == oldBrother) {
		    // secondNodeIndices.add(n);
		    // probabilities.add(backward);
		    // sum += backward;
		} else if ((i != j) && (tree.getNodeHeight(j) < heightIP)
			&& (heightIP < tree.getNodeHeight(jP))
			&& getNodeDistance(iP, jP) <= MAX_DISTANCE) {
		    secondNodeIndices.add(n);

		    pruneAndRegraft(tree, i, iP, j, jP);
		    double prob = Math.exp(calculateTreeLikelihood(prior,
			    likelihood, tree)
			    + offset);
		    probabilities.add(prob);
		    scores[n] = prob;
		    sum += prob;

		    pruneAndRegraft(tree, i, iP, oldBrother, oldGrandfather);
		}
	    }
	}

	if (sum <= 1E-100) {
	    // hack
	    // the proposals have such a small likelihood that they can be
	    // neglected
	    throw new OperatorFailedException(
		    "Couldn't find another proposal with a decent likelihood.");
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

	// int distance = getNodeDistance(iP, jP);
	// distances[distance]++;

	pruneAndRegraft(tree, i, iP, j, jP);

	// now simulate the backward move
	double sumBackward = 0.0;
	final NodeRef newBrother = j;
	final NodeRef newGrandfather = jP;
	for (int n = 0; n < nodeCount; n++) {
	    j = tree.getNode(n);
	    if (j != root) {
		jP = tree.getParent(j);

		if (j == newBrother) {
		    // secondNodeIndices.add(n);
		    // probabilities.add(backward);
		    // sum += backward;
		} else if ((i != j) && (tree.getNodeHeight(j) < heightIP)
			&& (heightIP < tree.getNodeHeight(jP))
			&& getNodeDistance(iP, jP) <= MAX_DISTANCE) {

		    if (scores[n] != Double.NEGATIVE_INFINITY) {
			sumBackward += scores[n];
		    } else {
			pruneAndRegraft(tree, i, iP, j, jP);
			double prob = Math.exp(calculateTreeLikelihood(prior,
				likelihood, tree)
				+ offset);
			sumBackward += prob;

			pruneAndRegraft(tree, i, iP, newBrother, newGrandfather);
		    }
		}
	    }
	}

	double forward = probabilities.get(index);

	// tree.pushTreeChangedEvent(jP);
	// tree.pushTreeChangedEvent(oldGrandfather);
	tree.pushTreeChangedEvent(i);

	double forwardProb = (forward / sum);
	double backwardProb = (backward / (sumBackward));
	double hastingsRatio = Math.log(backwardProb / forwardProb);

	return hastingsRatio;
    }

    private int getNodeDistance(NodeRef i, NodeRef j) {
	int count = 0;
	double heightI = tree.getNodeHeight(i);
	double heightJ = tree.getNodeHeight(j);

	while (i != j) {
	    count++;
	    if (heightI < heightJ) {
		i = tree.getParent(i);
		heightI = tree.getNodeHeight(i);
	    } else {
		j = tree.getParent(j);
		heightJ = tree.getNodeHeight(j);
	    }
	}
	return count;
    }

    public void printDistances() {
	System.out.println("Number of proposed trees in distances:");
	for (int i = 0; i < distances.length; i++) {
	    System.out.println(i + ")\t\t" + distances[i]);
	}
    }

    private double calculateTreeLikelihood(Prior prior, Likelihood likelihood,
	    TreeModel tree) {
	return evaluate(likelihood, prior);
    }

    private void pruneAndRegraft(TreeModel tree, NodeRef i, NodeRef iP,
	    NodeRef j, NodeRef jP) throws OperatorFailedException {
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
	    tree.endTreeEditUnsafe();
	} catch (MutableTree.InvalidTreeException ite) {
	    throw new OperatorFailedException(ite.toString());
	}

	tree.pushTreeChangedEvent(i);
    }

    /**
     * @param tree
     *            the tree
     * @param parent
     *            the parent
     * @param child
     *            the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

	if (tree.getChild(parent, 0) == child) {
	    return tree.getChild(parent, 1);
	} else {
	    return tree.getChild(parent, 0);
	}
    }

    /*
     * (non-Javadoc)
     * 
     * @see dr.evomodel.operators.SimpleGibbsOperator#getOperatorName()
     */
    @Override
    public String getOperatorName() {
	return GIBBS_PRUNE_AND_REGRAFT;
    }

    /*
     * (non-Javadoc)
     * 
     * @see dr.evomodel.operators.SimpleGibbsOperator#getStepCount()
     */
    @Override
    public int getStepCount() {
	return 0;
    }

    public static XMLObjectParser GIBBS_PRUNE_AND_REGRAFT_PARSER = new AbstractXMLObjectParser() {

	public String getParserName() {
	    return GIBBS_PRUNE_AND_REGRAFT;
	}

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	    TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
	    double weight = xo.getDoubleAttribute("weight");
	    boolean pruned = true;
	    if (xo.hasAttribute("pruned")) {
		pruned = xo.getBooleanAttribute("pruned");
	    }

	    return new GibbsPruneAndRegraft(treeModel, pruned, weight);
	}

	// ************************************************************************
	// AbstractXMLObjectParser implementation
	// ************************************************************************

	public String getParserDescription() {
	    return "This element represents a Gibbs sampler implemented through a prune and regraft operator. "
		    + "This operator prunes a random subtree and regrafts it below a node chosen by an importance distribution which is the proportion of the likelihoods of the proposals.";
	}

	public Class getReturnType() {
	    return GibbsPruneAndRegraft.class;
	}

	public XMLSyntaxRule[] getSyntaxRules() {
	    return rules;
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		AttributeRule.newDoubleRule("weight"),
		AttributeRule.newBooleanRule("pruned"),
		new ElementRule(TreeModel.class) };

    };
}
