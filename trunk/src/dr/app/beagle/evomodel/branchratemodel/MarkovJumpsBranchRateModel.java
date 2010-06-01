package dr.app.beagle.evomodel.branchratemodel;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * A branch rate model that takes a MarkovJumpsTreeLikelihood and gives a average of the rate
 * increments given by the product of rateParameter indicatorParameter weighted by the rewards
 * for each state.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MarkovJumpsBranchRateModel extends AbstractBranchRateModel {

    public MarkovJumpsBranchRateModel(String name, TreeModel treeModel, MarkovJumpsBeagleTreeLikelihood markovJumpsBeagleTreeLikelihood, Parameter rateParameter, Parameter indicatorParameter) {
        super(name);
        this.treeModel = treeModel;
        this.markovJumpsBeagleTreeLikelihood = markovJumpsBeagleTreeLikelihood;
        this.rateParameter = rateParameter;
        this.indicatorParameter = indicatorParameter;

        addModel(markovJumpsBeagleTreeLikelihood);
        addModel(treeModel);
        addVariable(rateParameter);
        addVariable(rateParameter);
    }

    public double getBranchRate(Tree tree, NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

        double[][] jumps = markovJumpsBeagleTreeLikelihood.getMarkovJumpsForNode(tree, node);

        double rate = 0.0;
        double total = 0.0;
        for (int i = 0; i < jumps.length; i++) {
            rate += rateParameter.getParameterValue(i) * indicatorParameter.getParameterValue(i) * jumps[i][0];
            total += indicatorParameter.getParameterValue(i) * jumps[i][0];
        }
        rate /= total;

        //System.out.println(rates[rateCategory] + "\t"  + rateCategory);
        return 1.0 + rate;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // this flags the entire tree for updating.
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
    }

    private final TreeModel treeModel;
    private final MarkovJumpsBeagleTreeLikelihood markovJumpsBeagleTreeLikelihood;
    private final Parameter rateParameter;
    private final Parameter indicatorParameter;

}
