package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.substmodel.DifferentialMassProvider;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
class BranchDifferentialMassProvider {

    private DifferentiableBranchRates indexHelper;
    private List<DifferentialMassProvider> differentialMassProviderList;

    BranchDifferentialMassProvider(DifferentiableBranchRates indexHelper,
                                   List<DifferentialMassProvider> differentialMassProviderList) {

        this.indexHelper = indexHelper;
        this.differentialMassProviderList = differentialMassProviderList;
    }

    double[] getDifferentialMassMatrixForBranch(NodeRef node, double time) {
        int index = indexHelper == null ? 0 : indexHelper.getParameterIndexFromNode(node);
        return differentialMassProviderList.get(index).getDifferentialMassMatrix(time);
    }
}
