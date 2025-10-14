package dr.evomodel.branchratemodel.latentStateBranchRate;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.SericolaLatentStateBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.Bounds;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import  dr.inference.model.Parameter;
/**
 * This statistic reports the length a tree spent either replicating or in the latent state.
 * It is a parameter so it can be used in distributions like the gamma distribution.
 * @author JT McCrone
 * 
 */


public class LatentStateBranchRateLengthStatistic extends Parameter.Abstract implements ModelListener {
    public LatentStateBranchRateLengthStatistic( SericolaLatentStateBranchRateModel model, Tree tree, STATE state) {
      
        this.model = model;
       
        this.state = state;
        this.tree= tree;
         model.addModelListener(this);
        if(tree instanceof TreeModel){
            ((TreeModel) tree).addModelListener(this);
        }
    }


    public int getDimension() {
        return 1;
    }

    /**
     * @return the total length of all the branches in the tree
     */
    public double getParameterValue(int dim) {
        double length =0.0;
        for(int i=0; i< tree.getNodeCount(); i++) {
            if(!tree.isRoot(tree.getNode(i))){
                double branchLength =tree.getBranchLength(tree.getNode(i));
                double latentProportion = model.getLatentProportion(tree, tree.getNode(i));
                if(state == STATE.REPLICATING) {
                    length += branchLength * (1-latentProportion);
                } else if(state == STATE.LATENT) {
                    length += branchLength * latentProportion;
                }
            }
        }
        return length;
    }

    @Override
    public boolean isImmutable() {
        return true;
    }

     protected void storeValues() {
      //do nothing
    }

    protected void restoreValues() {
      //do nothing
    }

    protected void acceptValues() {
      //do nothing
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }


    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        return getId();
    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }
    @Override
     public void modelChangedEvent(Model model, Object object, int index) {
       fireParameterChangedEvent();
     }


     @Override
     public void modelRestored(Model model) {
            // do nothing
     }

    private SericolaLatentStateBranchRateModel model = null;
    private STATE state;
    private Tree tree;
    public static String LATENT_STATE_BRANCH_RATE_LENGTH_STATISTIC = "latentStateBranchRateLengthStatistic";
    public enum STATE {REPLICATING, LATENT}
     private Bounds bounds = null;

}
