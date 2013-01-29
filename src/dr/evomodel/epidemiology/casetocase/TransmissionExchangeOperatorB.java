package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

/**
 * Implements branch exchange operations that also exchange entire subtrees of the transmission tree. As this already
 * severely restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange. Needs a
 * parser.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorB extends AbstractTreeOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;

    public TransmissionExchangeOperatorB(CaseToCaseTransmissionLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();

        final int tipCount = tree.getExternalNodeCount();

        assert tree.getExternalNodeCount() == tipCount :
                "Lost some tips";

        return 0;
    }

    public void exchange() throws OperatorFailedException{
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;
        NodeRef iP = tree.getParent(i);
        boolean paintingsDoNotMatch = false;

        while(root == i || !paintingsDoNotMatch) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
            iP = tree.getParent(i);
            if(root!=i){
                paintingsDoNotMatch = branchMap[i.getNumber()]!=branchMap[iP.getNumber()];
            }
        }

        NodeRef j = root;
        NodeRef jP = tree.getParent(j);
        paintingsDoNotMatch = false;

        while(root == j || i==j || iP==jP || !paintingsDoNotMatch) {
            j = tree.getNode(MathUtils.nextInt(nodeCount));
            jP = tree.getParent(j);
            if(root!=j){
                paintingsDoNotMatch = branchMap[j.getNumber()]!=branchMap[jP.getNumber()];
            }
        }

/*
        Intuitively it would seem this is a lot more likely to succeed than operator A.
*/

        if((i != jP) && (j != iP)
                && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                && (tree.getNodeHeight(i) < tree.getNodeHeight(jP)) ) {
            exchangeNodes(tree, i, j, iP, jP);
            return;
        }

        // What's happened to the node numbers? If the node has its old number in its new position, you don't need
        // to even update the branch map!

        c2cLikelihood.makeDirty();

        throw new OperatorFailedException("Couldn't find valid move on this tree!");

    }



    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return null;
    }
}
