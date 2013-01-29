package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Implements branch exchange operations that leave the transmission network unchanged. As this already severely
 * restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange. Needs a parser.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorA extends AbstractTreeOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;

    public TransmissionExchangeOperatorA(CaseToCaseTransmissionLikelihood c2cLikelihood, double weight) {
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

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;
        NodeRef iP = tree.getParent(i);
        HashSet<Integer> possibleSwaps = c2cLikelihood.samePainting(iP,false);
        int noPossibleSwaps = 0;

        while(root == i || noPossibleSwaps == 1) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
            iP = tree.getParent(i);
            possibleSwaps = c2cLikelihood.samePainting(iP,false);
            noPossibleSwaps = possibleSwaps.size();
        }

        NodeRef jP=iP;
        ArrayList<Integer> swapsAsList = new ArrayList<Integer>(possibleSwaps);

        while(jP==iP){
            Collections.shuffle(swapsAsList);
            jP = tree.getNode(swapsAsList.get(0));
        }

        NodeRef j = tree.getChild(jP, MathUtils.nextInt(1));

/*
        I tend to think that this may fail quite a lot of the time due to lack of candidates... a version that does
        actually adjust heights might be necessary in the long run. Narrow exchange might be much more likely to
        actually succeed in changing the tree if the paintings allow the tree to be changed in that way; might
        have to investigate which problem is more serious.
*/

        if( (i != jP) && (j != iP)
                && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                && (tree.getNodeHeight(i) < tree.getNodeHeight(jP)) ) {
            exchangeNodes(tree, i, j, iP, jP);
            return;
        }

        // What's happened to the node numbers?



        c2cLikelihood.makeDirty();

        throw new OperatorFailedException("Couldn't find valid wide move on this tree!");

    }



    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return null;
    }

}
