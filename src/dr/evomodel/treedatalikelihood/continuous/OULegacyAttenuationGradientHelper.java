package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

public final class OULegacyAttenuationGradientHelper {

    private OULegacyAttenuationGradientHelper() {
    }

    public static DenseMatrix64F getVarianceWrtAttenuationDiagonal(final OUDiffusionModelDelegate delegate,
                                                                   final NodeRef node,
                                                                   final ContinuousDiffusionIntegrator cdi,
                                                                   final BranchSufficientStatistics statistics,
                                                                   final DenseMatrix64F gradient) {
        return OULegacyActualizationGradientHelper.getGradientVarianceWrtAttenuationDiagonal(
                delegate.getElasticModel(),
                delegate.getDimTrait(),
                cdi,
                delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                delegate.getEigenBufferOffsetIndex(0),
                statistics,
                node.getNumber(),
                gradient);
    }

    public static DenseMatrix64F getDisplacementWrtAttenuationDiagonal(final OUDiffusionModelDelegate delegate,
                                                                       final NodeRef node,
                                                                       final ContinuousDiffusionIntegrator cdi,
                                                                       final BranchSufficientStatistics statistics,
                                                                       final DenseMatrix64F gradient) {
        return OULegacyActualizationGradientHelper.getGradientDisplacementWrtAttenuationDiagonal(
                delegate.getDimTrait(),
                cdi,
                delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                statistics,
                delegate.getBranchDriftRate(node),
                gradient);
    }

    public static DenseMatrix64F unsupportedVarianceWrtAttenuation(final String kind) {
        throw new UnsupportedOperationException(
                "Legacy attenuation gradient is not implemented for " + kind +
                        " OU actualization; use the canonical branch path for that family.");
    }

    public static DenseMatrix64F unsupportedDisplacementWrtAttenuation(final String kind) {
        throw new UnsupportedOperationException(
                "Legacy displacement-vs-attenuation gradient is not implemented for " + kind +
                        " OU actualization; use the canonical branch path for that family.");
    }
}
