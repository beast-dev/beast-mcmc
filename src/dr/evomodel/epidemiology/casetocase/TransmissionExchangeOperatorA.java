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
 * restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.
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
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;

        while( root == i ) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        }

        NodeRef iP = tree.getParent(i);

        HashSet<Integer> possibleSwaps = c2cLikelihood.samePainting(iP,false);

        NodeRef jP=iP;
        ArrayList<Integer> asList = new ArrayList<Integer>(possibleSwaps);

        while(jP==iP){
            Collections.shuffle(asList);
            jP = tree.getNode(asList.get(0));
        }

        NodeRef j = i;

        while(j==i || j==root){

        }


    }



    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

}
