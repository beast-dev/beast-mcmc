package dr.app.beagle.evomodel.operators;

import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractLikelihoodCore;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TipStateSwapOperator extends SimpleMCMCOperator {

    public static final String TIP_STATE_OPERATOR = dr.evomodel.operators.TipStateSwapOperator.TIP_STATE_OPERATOR;

    public TipStateSwapOperator(AncestralStateBeagleTreeLikelihood treeLikelihood, double weight) {
        this.treeLikelihood = treeLikelihood;
        setWeight(weight);
        int patternCount = treeLikelihood.getPatternCount();
        states1 = new int[patternCount];
        states2 = new int[patternCount];
    }

    private final int[] states1;
    private final int[] states2;
    private int index1;
    private int index2;

    public double doOperation() throws OperatorFailedException {

        int tipCount = treeLikelihood.getTreeModel().getExternalNodeCount();

        // Choose two tips to swap
        index1 = MathUtils.nextInt(tipCount);
        index2 = index1;
        while (index2 == index1)
            index2 = MathUtils.nextInt(tipCount);

        swap(index1, index2);

        treeLikelihood.makeDirty();

        return 0;
    }

    private void swap(int i, int j) {

        treeLikelihood.getStates(i, states1);
        treeLikelihood.getStates(j, states2);

        treeLikelihood.setStates(j, states1);
        treeLikelihood.setStates(i, states2);                
    }

    public void reject() {
        super.reject();
        // There is currently no restore functions for tip states, so manually adjust state
        swap(index1, index2);
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String getOperatorName() {
        return TIP_STATE_OPERATOR;
    }

    private final AncestralStateBeagleTreeLikelihood treeLikelihood;
}