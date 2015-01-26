package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.FixedDriftModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Created by mandevgill on 12/22/14.
 */
public class FixedDriftModel extends AbstractBranchRateModel {

    public FixedDriftModel(Parameter rateOne,
                           Parameter rateTwo,
                           Parameter remainingRates,
                           String rateOneID,
                           String rateTwoID) {

        super(FixedDriftModelParser.FIXED_DRIFT);


        //   addModel(treeModel);
        //  this.treeModel = treeModel;

        this.rateOne = rateOne;
        addVariable(rateOne);

        this.rateTwo = rateTwo;
        addVariable(rateTwo);

        this.remainingRates = remainingRates;
        addVariable(remainingRates);

        this.rateOneID = rateOneID;
        this.rateTwoID = rateTwoID;
    }


    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        //  calculateBranchRates(treeModel);
        // recalculateScaleFactor();
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        if (tree.isExternal(node)) {
            if (tree.getTaxonId(node.getNumber()).matches(rateOneID)) {
                //  System.err.println("node one id " + tree.getTaxonId(node.getNumber()));
                return rateOne.getParameterValue(0);
            } else if (tree.getTaxonId(node.getNumber()).matches(rateTwoID)) {
                // System.err.println("node two id " + tree.getTaxonId(node.getNumber()));
                return rateTwo.getParameterValue(0);
            } else {
                return remainingRates.getParameterValue(0);
            }
        } else {
            return remainingRates.getParameterValue(0);
        }
    }


    // the tree model
    private TreeModel treeModel;

    // the unscaled rates of each branch, taking into account the indicators
    private double[] branchRates;


    private Parameter remainingRates;
    private Parameter rateOne;
    private Parameter rateTwo;
    private String rateOneID;
    private String rateTwoID;

}
