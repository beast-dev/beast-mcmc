/**
 *
 */
package dr.evomodel.operators;

import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMetropolizedGibbsOperator;
import dr.inference.prior.Prior;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Sebastian Hoehna
 *
 */
public class GibbsSubtreeSwap extends SimpleMetropolizedGibbsOperator {

    private int MAX_DISTANCE = 10;

    private final TreeModel tree;

    private final int[] distances;

    private boolean pruned = true;

    /**
     *
     */
    public GibbsSubtreeSwap(TreeModel tree, boolean pruned, double weight) {
        this.tree = tree;
        this.pruned = pruned;
        setWeight(weight);
        MAX_DISTANCE = tree.getNodeCount() / 10;
        MAX_DISTANCE = 4;
        distances = new int[tree.getNodeCount() / 2];
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

    public double doOperation(Prior prior, Likelihood likelihood)
            throws OperatorFailedException {

        if( pruned ) {
            return prunedWide(prior, likelihood);
        } else {
            return wide(prior, likelihood);
        }
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     *
     * @throws InvalidTreeException
     */
    public double wide(Prior prior, Likelihood likelihood)
            throws OperatorFailedException {

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i;
        int indexI;

        do {
            indexI = MathUtils.nextInt(nodeCount);
            i = tree.getNode(indexI);
        } while( root == i
                || (tree.getParent(i) == root && tree.getNodeHeight(i) > tree
                .getNodeHeight(getOtherChild(tree, tree.getParent(i), i))) );

        List<Integer> secondNodeIndices = new ArrayList<Integer>();
        List<Double> probabilities = new ArrayList<Double>();
        NodeRef j, jP;
        NodeRef iP = tree.getParent(i);
        double heightIP = tree.getNodeHeight(iP);
        double heightI = tree.getNodeHeight(i);
        double sum = 0.0;
        double backward = calculateTreeLikelihood(prior, likelihood, tree);
        int offset = (int) -backward;
        backward = Math.exp(backward + offset);
        for(int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if( j != root ) {
                jP = tree.getParent(j);

                if( (iP != jP) && (tree.getNodeHeight(j) < heightIP)
                        && (heightI < tree.getNodeHeight(jP)) ) {
                    secondNodeIndices.add(n);

                    swap(tree, i, j, iP, jP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    probabilities.add(prob);
                    swap(tree, i, j, jP, iP);
                    sum += prob;

                }
            }
        }

        if( sum <= 1E-100 ) {
            // hack
            // the proposals have such a small likelihood that they can be
            // neglected
            throw new OperatorFailedException(
                    "Couldn't find another proposal with a decent likelihood.");
        }

        double ran = Math.random() * sum;
        int index = 0;
        while( ran > 0.0 ) {
            ran -= probabilities.get(index);
            index++;
        }
        index--;

        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);
        double heightJP = tree.getNodeHeight(jP);
        double heightJ = tree.getNodeHeight(j);

        // int distance = getNodeDistance(i, j);
        // distances[distance]++;

        // *******************************************
        // assuming we would have chosen j first
        double sumForward2 = 0.0;
        NodeRef k, kP;
        for(int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if( k != root ) {
                kP = tree.getParent(k);

                if( (jP != kP) && (tree.getNodeHeight(k) < heightJP)
                        && (heightJ < tree.getNodeHeight(kP)) ) {

                    swap(tree, j, k, jP, kP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumForward2 += prob;
                    swap(tree, j, k, kP, jP);
                }
            }
        }

        swap(tree, i, j, iP, jP);
        double forward = probabilities.get(index);

        iP = jP;
        heightIP = heightJP;
        double sumBackward = 0.0;
        for(int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if( j != root ) {
                jP = tree.getParent(j);

                if( (iP != jP) && (tree.getNodeHeight(j) < heightIP)
                        && (heightI < tree.getNodeHeight(jP)) ) {

                    swap(tree, i, j, iP, jP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumBackward += prob;
                    swap(tree, i, j, jP, iP);

                }
            }
        }

        // *******************************************
        // assuming we would have chosen j first
        double sumBackward2 = 0.0;
        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);
        heightJP = tree.getNodeHeight(jP);
        heightJ = tree.getNodeHeight(j);
        for(int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if( k != root ) {
                kP = tree.getParent(k);

                if( (jP != kP) && (tree.getNodeHeight(k) < heightJP)
                        && (heightJ < tree.getNodeHeight(kP)) ) {

                    swap(tree, j, k, jP, kP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumBackward2 += prob;
                    swap(tree, j, k, kP, jP);
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

    /**
     * WARNING: Assumes strictly bifurcating tree.
     *
     * @throws InvalidTreeException
     */
    public double prunedWide(Prior prior, Likelihood likelihood)
            throws OperatorFailedException {

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i;
        int indexI;

        do {
            indexI = MathUtils.nextInt(nodeCount);
            i = tree.getNode(indexI);
        } while( root == i
                || (tree.getParent(i) == root && tree.getNodeHeight(i) > tree
                .getNodeHeight(getOtherChild(tree, tree.getParent(i), i))) );

        List<Integer> secondNodeIndices = new ArrayList<Integer>();
        List<Double> probabilities = new ArrayList<Double>();
        NodeRef j, jP;
        NodeRef iP = tree.getParent(i);
        double heightIP = tree.getNodeHeight(iP);
        double heightI = tree.getNodeHeight(i);
        double sum = 0.0;
        double backward = calculateTreeLikelihood(prior, likelihood, tree);
        int offset = (int) -backward;
        backward = Math.exp(backward + offset);
        for(int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if( j != root ) {
                jP = tree.getParent(j);

                if( (iP != jP) && (tree.getNodeHeight(j) < heightIP)
                        && (heightI < tree.getNodeHeight(jP))
                        && getNodeDistance(iP, jP) <= MAX_DISTANCE ) {
                    secondNodeIndices.add(n);

                    swap(tree, i, j, iP, jP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    probabilities.add(prob);
                    swap(tree, i, j, jP, iP);
                    sum += prob;

                }
            }
        }

        if( sum <= 1E-100 ) {
            // hack
            // the proposals have such a small likelihood that they can be
            // neglected
            throw new OperatorFailedException(
                    "Couldn't find another proposal with a decent likelihood.");
        }

        double ran = Math.random() * sum;
        int index = 0;
        while( ran > 0.0 ) {
            ran -= probabilities.get(index);
            index++;
        }
        index--;

        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);
        double heightJP = tree.getNodeHeight(jP);
        double heightJ = tree.getNodeHeight(j);

        // *******************************************
        // assuming we would have chosen j first
        double sumForward2 = 0.0;
        NodeRef k, kP;
        for(int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if( k != root ) {
                kP = tree.getParent(k);

                if( (jP != kP) && (tree.getNodeHeight(k) < heightJP)
                        && (heightJ < tree.getNodeHeight(kP))
                        && getNodeDistance(kP, jP) <= MAX_DISTANCE ) {

                    swap(tree, j, k, jP, kP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumForward2 += prob;
                    swap(tree, j, k, kP, jP);
                }
            }
        }

        swap(tree, i, j, iP, jP);
        double forward = probabilities.get(index);

        iP = jP;
        heightIP = heightJP;
        double sumBackward = 0.0;
        for(int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if( j != root ) {
                jP = tree.getParent(j);

                if( (iP != jP) && (tree.getNodeHeight(j) < heightIP)
                        && (heightI < tree.getNodeHeight(jP))
                        && getNodeDistance(iP, jP) <= MAX_DISTANCE ) {

                    swap(tree, i, j, iP, jP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumBackward += prob;
                    swap(tree, i, j, jP, iP);

                }
            }
        }

        // *******************************************
        // assuming we would have chosen j first
        double sumBackward2 = 0.0;
        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);
        heightJP = tree.getNodeHeight(jP);
        heightJ = tree.getNodeHeight(j);
        for(int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if( k != root ) {
                kP = tree.getParent(k);

                if( (jP != kP) && (tree.getNodeHeight(k) < heightJP)
                        && (heightJ < tree.getNodeHeight(kP))
                        && getNodeDistance(kP, jP) <= MAX_DISTANCE ) {

                    swap(tree, j, k, jP, kP);
                    double prob = Math.exp(calculateTreeLikelihood(prior,
                            likelihood, tree)
                            + offset);
                    sumBackward2 += prob;
                    swap(tree, j, k, kP, jP);
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

    private double calculateTreeLikelihood(Prior prior, Likelihood likelihood,
                                           TreeModel tree) {
        return evaluate(likelihood, prior, 1.0);
        // return 0.0;
    }

    /**
     * @param tree   the tree
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

        if( tree.getChild(parent, 0) == child ) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
    }

    /* exchange subtrees whose root are i and j */
    private TreeModel swap(TreeModel tree, NodeRef i, NodeRef j, NodeRef iP,
                           NodeRef jP) throws OperatorFailedException {

        tree.beginTreeEdit();
        tree.removeChild(iP, i);
        tree.removeChild(jP, j);
        tree.addChild(jP, i);
        tree.addChild(iP, j);

        tree.endTreeEdit();

        return tree;
    }

    private int getNodeDistance(NodeRef i, NodeRef j) {
        int count = 0;
        double heightI = tree.getNodeHeight(i);
        double heightJ = tree.getNodeHeight(j);

        while( i != j ) {
            count++;
            if( heightI < heightJ ) {
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
        for(int i = 0; i < distances.length; i++) {
            System.out.println(i + ")\t\t" + distances[i]);
        }
    }

    public String getOperatorName() {
        return "Gibbs Subtree Exchange";
    }

}
