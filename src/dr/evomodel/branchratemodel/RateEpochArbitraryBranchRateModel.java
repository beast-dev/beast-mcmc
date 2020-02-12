package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.branchratemodel.RateEpochArbitraryBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class RateEpochArbitraryBranchRateModel extends AbstractBranchRateModel {

    private final RateEpochBranchRateModel epochBranchRateModel;
    private final ArbitraryBranchRates arbitraryBranchRates;

    public RateEpochArbitraryBranchRateModel(RateEpochBranchRateModel epochBranchRateModel,
                                             ArbitraryBranchRates arbitraryBranchRates) {
        super(RateEpochArbitraryBranchRateModelParser.RATE_EPOCH_ARBITRARY_BRANCH_RATES);

        this.epochBranchRateModel = epochBranchRateModel;
        this.arbitraryBranchRates = arbitraryBranchRates;

        addModel(epochBranchRateModel);
        addModel(arbitraryBranchRates);
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        return epochBranchRateModel.getBranchRate(tree, node) * arbitraryBranchRates.getBranchRate(tree, node);
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
