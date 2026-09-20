package dr.evomodel.epidemiology;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;

public class TwoPathogenModelLikelihood extends AbstractModelLikelihood {

    private final TwoPathogenModel twoPathogenModel;
    private boolean likelihoodKnown = false;
    private final CompartmentalModelSimulator simulator;
    private final TreeModel treeModelOne;
    private final TreeModel treeModelTwo;
    private final double mostRecentSamplingDateOne;
    private final double mostRecentSamplingDateTwo;
    private final int numGridPoints;
    private final double intervalWidth;

    public TwoPathogenModelLikelihood(TwoPathogenModel twoPathogenModel,
                                      CompartmentalModelSimulator simulator,
                                      TreeModel treeModelOne,
                                      TreeModel treeModelTwo) {

        super("TwoPathogenModelLikelihood");

        this.twoPathogenModel = twoPathogenModel;
        this.simulator = simulator;
        this.treeModelOne = treeModelOne;
        this.treeModelTwo = treeModelTwo;
        this.mostRecentSamplingDateOne = twoPathogenModel.mostRecentSamplingDateOne;
        this.mostRecentSamplingDateTwo = twoPathogenModel.mostRecentSamplingDateTwo;
        this.numGridPoints = twoPathogenModel.numGridPoints;
        this.intervalWidth = twoPathogenModel.cutOff / numGridPoints;
        addModel(twoPathogenModel);
        addModel(treeModelOne);
        addModel(treeModelTwo);
        // make sure initial lineage counts are provided to simulator
        updateLineageCounts();

        // run stochastic simulator until we get a valid initial trajectory
        int maxAttempts = 100000;
        int attempts = 0;
        do{
            simulator.simulateTrajectory();
            attempts++;
            if(attempts % 1000 == 0){
                System.out.println("attempting to find valid initial trajectory: attempt " + attempts);
            }
            if(attempts >= maxAttempts){
                throw new RuntimeException("Could not find valid initial trajectory of compartment " +
                        "counts after " + maxAttempts + " attempts. Check that model parameters are consistent with trees.");
            }

        }while(simulator.isLineageConstraintViolated());
        System.out.println("Found valid initial trajectory after " + attempts + " attempts.");
    }

    private void updateLineageCounts(){
        int[][] lineageCounts = new int[2][numGridPoints];
        lineageCounts[0] = computeLineageCounts(treeModelOne, mostRecentSamplingDateOne);
        lineageCounts[1] = computeLineageCounts(treeModelTwo, mostRecentSamplingDateTwo);
        simulator.setLineageCounts(lineageCounts);

        System.out.println("Lineage counts:");
        for (int k = 0; k < numGridPoints; k++) {
            System.out.println("Grid point " + k + ": pathogen1=" +
                    lineageCounts[0][k] + " pathogen2=" + lineageCounts[1][k]);
        }
    }

    private int[] computeLineageCounts(TreeModel treeModel, double mostRecentSamplingDate){
        int[] lineageCounts = new int[numGridPoints];
        double mostRecentDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        double cutOff = twoPathogenModel.cutOff;

        for(int k = 0; k < numGridPoints; k++){
            // backward time from this tree's most recent sampling date
            double backwardTime = cutOff - k*intervalWidth - (mostRecentDate-mostRecentSamplingDate);

            // backwardTime being negative means that the grid point is more recent than the tree's
            // most recent sampling date, in which case lineageCounts should be set to 0
            if(backwardTime < 0){
                lineageCounts[k] = 0;
            }else{
                lineageCounts[k] = countLineagesAtTime(treeModel, backwardTime);
            }
        }
        return lineageCounts;
    }

    // Count the number of lineages in a tree at the given backard time
    // check if there is a more efficient way to do this
    private int countLineagesAtTime(TreeModel treeModel, double backwardTime){
        int count = 0;
        for (int i = 0; i < treeModel.getNodeCount(); i++){
            NodeRef node = treeModel.getNode(i);
            if(!treeModel.isRoot(node)){
                double nodeHeight = treeModel.getNodeHeight(node);
                double parentHeight = treeModel.getNodeHeight(treeModel.getParent(node));
                if(backwardTime >= nodeHeight && backwardTime < parentHeight){
                    count++;
                }
            }
        }
        return count;
    }

    public double getLogLikelihood() {
        // return -infinity if lineage count constraint is violated
        if (simulator.isLineageConstraintViolated()){
            return Double.NEGATIVE_INFINITY;
        }
        // otherwise, always return 0.0 since trajectory prior and proposal cancel in MH ratio
        return 0.0;
    }

    public void makeDirty(){
        // do nothing
    }

    protected void handleModelChangedEvent(Model model, Object object, int index){
        if(model == treeModelOne || model == treeModelTwo){
            updateLineageCounts();
        }
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState(){
        // do nothing
    }

    protected void restoreState(){
        simulator.resetLineageConstraintViolated();
    }

    protected void acceptState(){
        // do nothing
    }

    public Model getModel(){
        return this;
    }

}
