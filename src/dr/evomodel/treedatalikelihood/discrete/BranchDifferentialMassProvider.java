package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.tree.TreeParameterModel;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
class BranchDifferentialMassProvider {

    private TreeParameterModel indexHelper;
    private List<DifferentialMassProvider> differentialMassProviderList;

    BranchDifferentialMassProvider(TreeParameterModel indexHelper,
                                   List<DifferentialMassProvider> differentialMassProviderList) {

        this.indexHelper = indexHelper;
        this.differentialMassProviderList = differentialMassProviderList;

    }

    double[] getDifferentialMassMatrixForBranch(NodeRef node, double time) {
        return differentialMassProviderList.get(indexHelper.getParameterIndexFromNodeNumber(node.getNumber())).getDifferentialMassMatrix(time);
    }
}
