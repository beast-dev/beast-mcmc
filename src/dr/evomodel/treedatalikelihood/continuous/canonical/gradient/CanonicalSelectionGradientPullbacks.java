package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;

final class CanonicalSelectionGradientPullbacks {

    private CanonicalSelectionGradientPullbacks() { }

    static CanonicalSelectionGradientPullback create(final OUProcessModel processModel,
                                                     final int dimension,
                                                     final double[] gradA,
                                                     final BranchGradientWorkspace workspace) {
        final SelectionMatrixParameterization parameterization =
                processModel.getSelectionMatrixParameterization();
        if (parameterization instanceof SpecializedCanonicalSelectionParameterization) {
            return new SpecializedCanonicalSelectionGradientPullback(
                    (SpecializedCanonicalSelectionParameterization) parameterization,
                    dimension,
                    gradA,
                    workspace);
        }
        return new DenseCanonicalSelectionGradientPullback(dimension);
    }
}
