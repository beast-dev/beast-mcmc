package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

import java.util.ArrayList;

/**
 * An operator that wraps a (phylogenetic) tree operator and adjusts the transmission tree accordingly. Only works
 * for the simpler (TT=TMRCA) version of the transmission tree. TT moves are deterministic, so the Hastings ratios
 * reported by the other operators are not affected.
 *
 * @author Matthew Hall
 */

public class TransmissionTreeOperator extends AbstractCoercableOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;
    private final AbstractTreeOperator innerOperator;

    public TransmissionTreeOperator(CaseToCaseTransmissionLikelihood c2cLikelihood, AbstractTreeOperator operator,
                                    CoercionMode mode) {
        super(mode);
        this.c2cLikelihood = c2cLikelihood;
        this.innerOperator = operator;
    }

    public TransmissionTreeOperator(CaseToCaseTransmissionLikelihood c2cLikelihood, AbstractTreeOperator operator) {
        this(c2cLikelihood,operator,CoercionMode.COERCION_OFF);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();
        int[] oldParents = getParentsArray(tree);
        double hr = innerOperator.doOperation();
        int[] newParents = getParentsArray(tree);
        ArrayList<Integer> changedNodes = new ArrayList<Integer>();
        for(int i=0; i<tree.getNodeCount(); i++){
            if(oldParents[i]!=newParents[i]){
                changedNodes.add(i);
            }
        }
        if(changedNodes.size()!=0){
            if(changedNodes.size()==2){
                //this is a node swap operator
            } else if(changedNodes.size()==3){
                //this is a node transplantation operator
            } else {
                //I don't know what this is
            }
        }
        return hr;
    }

    private int[] getParentsArray(TreeModel tree){
        int[] out = new int[tree.getNodeCount()];
        for(int i=0; i<tree.getNodeCount(); i++){
            if(tree.getNode(i)==tree.getRoot()){
                out[i]=-1;
            } else {
                out[i]=tree.getParent(tree.getNode(i)).getNumber();
            }
        }
        return out;
    }


    @Override
    public double getCoercableParameter() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setCoercableParameter(double value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getRawParameter() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
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
