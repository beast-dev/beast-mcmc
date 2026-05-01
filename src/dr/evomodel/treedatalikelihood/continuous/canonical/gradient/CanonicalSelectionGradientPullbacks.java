package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;

final class CanonicalSelectionGradientPullbacks {

    private CanonicalSelectionGradientPullbacks() { }

    static CanonicalSelectionGradientPullback create(final OUProcessModel processModel,
                                                     final int dimension,
                                                     final double[] gradA,
                                                     final BranchGradientWorkspace workspace) {
        if (processModel.getSelectionMatrixParameterization()
                instanceof SpecializedCanonicalSelectionParameterization) {
            return new SpecializedCanonicalSelectionGradientPullback(
                    (SpecializedCanonicalSelectionParameterization)
                            processModel.getSelectionMatrixParameterization(),
                    dimension,
                    gradA,
                    workspace);
        }
        return new DenseCanonicalSelectionGradientPullback(dimension);
    }
}
