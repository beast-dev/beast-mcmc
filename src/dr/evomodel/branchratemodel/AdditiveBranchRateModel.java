package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.branchratemodel.AdditiveBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class AdditiveBranchRateModel extends AbstractBranchRateModel {

    private final List<AbstractBranchRateModel> branchRateModels;
    private final int dim;

    public AdditiveBranchRateModel(List<AbstractBranchRateModel> branchRateModels) {
        super(AdditiveBranchRateModelParser.ADDITIVE_BRANCH_RATES);

        this.branchRateModels = branchRateModels;
        this.dim = branchRateModels.size();

        for (AbstractBranchRateModel branchRateModel : branchRateModels) {
            addModel(branchRateModel);
        }
    }
    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        double rate = branchRateModels.get(0).getBranchRate(tree, node);
        for (int i = 1; i < dim; i++) {
            rate += branchRateModels.get(i).getBranchRate(tree, node);
        }
        return rate;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
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

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }
}
