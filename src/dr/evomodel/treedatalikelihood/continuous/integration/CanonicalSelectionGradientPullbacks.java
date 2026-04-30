package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;

final class CanonicalSelectionGradientPullbacks {

    private CanonicalSelectionGradientPullbacks() { }

    static CanonicalSelectionGradientPullback create(final OUProcessModel processModel,
                                                     final int dimension,
                                                     final double[] gradA,
                                                     final BranchGradientWorkspace workspace) {
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            return new OrthogonalBlockCanonicalSelectionGradientPullback(
                    (OrthogonalBlockCanonicalParameterization)
                            processModel.getSelectionMatrixParameterization(),
                    dimension,
                    gradA,
                    workspace);
        }
        return new DenseCanonicalSelectionGradientPullback(dimension);
    }
}
