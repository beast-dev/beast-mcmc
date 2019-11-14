package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.tree.TreeParameterModel;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
class BranchDifferentialMassProvider {

    private IndexHelper indexHelper;
    private List<DifferentialMassProvider> differentialMassProviderList;

    BranchDifferentialMassProvider(TreeParameterModel indexHelper,
                                   List<DifferentialMassProvider> differentialMassProviderList) {

        this.indexHelper = factory(indexHelper);
        this.differentialMassProviderList = differentialMassProviderList;

    }

    double[] getDifferentialMassMatrixForBranch(NodeRef node, double time) {
        return differentialMassProviderList.get(indexHelper.getParameterIndexFromNodeNumber(node.getNumber())).getDifferentialMassMatrix(time);
    }

    private IndexHelper factory(TreeParameterModel indexHelper) {
        if (indexHelper == null) {
            return new IndexHelper.HomogeneousIndexHelper();
        } else {
            return new IndexHelper.TreeParameterIndexHelper(indexHelper);
        }
    }

    private interface IndexHelper {
        int getParameterIndexFromNodeNumber(int nodeNumber);

        class TreeParameterIndexHelper implements IndexHelper {

            private TreeParameterModel helper;

            TreeParameterIndexHelper(TreeParameterModel helper) {
                this.helper = helper;
            }

            @Override
            public int getParameterIndexFromNodeNumber(int nodeNumber) {
                return helper.getParameterIndexFromNodeNumber(nodeNumber);
            }
        }

        class HomogeneousIndexHelper implements IndexHelper {

            @Override
            public int getParameterIndexFromNodeNumber(int nodeNumber) {
                return 0;
            }
        }
    }
}
