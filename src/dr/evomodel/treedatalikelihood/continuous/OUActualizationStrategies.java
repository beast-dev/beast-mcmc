package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.unwrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

final class OUActualizationStrategyFactory {

    private OUActualizationStrategyFactory() {
    }

    static OUActualizationStrategy create(final OUDelegateContext context) {
        if (context.hasDiagonalActualization()) {
            return new OUActualizationStrategyDiagonal(context);
        }
        if (context.hasBlockDiagActualization()) {
            return new OUActualizationStrategyBlockDecomposition(context);
        }
        return new OUActualizationStrategyGeneral(context);
    }
}

abstract class AbstractOUActualizationStrategy implements OUActualizationStrategy {

    protected final OUDelegateContext context;

    AbstractOUActualizationStrategy(final OUDelegateContext context) {
        this.context = context;
    }

    @Override
    public DenseMatrix64F gradientVarianceWrtVariance(final NodeRef node,
                                                      final ContinuousDiffusionIntegrator cdi,
                                                      final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                      final DenseMatrix64F gradient) {
        throw new UnsupportedOperationException("Gradient variance path not yet ported in local CDI prototype");
    }

    @Override
    public DenseMatrix64F gradientVarianceWrtAttenuation(final NodeRef node,
                                                         final ContinuousDiffusionIntegrator cdi,
                                                         final BranchSufficientStatistics statistics,
                                                         final DenseMatrix64F dSigma) {
        throw new UnsupportedOperationException("Selection variance adjoint not yet ported in local CDI prototype");
    }

    @Override
    public DenseMatrix64F gradientDisplacementWrtDrift(final NodeRef node,
                                                       final ContinuousDiffusionIntegrator cdi,
                                                       final DenseMatrix64F gradient) {
        throw new UnsupportedOperationException("Drift displacement adjoint not yet ported in local CDI prototype");
    }

    @Override
    public DenseMatrix64F gradientDisplacementWrtAttenuation(final NodeRef node,
                                                             final ContinuousDiffusionIntegrator cdi,
                                                             final BranchSufficientStatistics statistics,
                                                             final DenseMatrix64F gradient) {
        throw new UnsupportedOperationException("Selection displacement adjoint not yet ported in local CDI prototype");
    }

    @Override
    public void meanTipVariances(final double priorSampleSize,
                                 final double[] treeLengths,
                                 final DenseMatrix64F traitVariance,
                                 final DenseMatrix64F varSum) {
        throw new UnsupportedOperationException("meanTipVariances not yet ported in local CDI prototype");
    }

    @Override
    public double[] rootGradient(final int index,
                                 final ContinuousDiffusionIntegrator cdi,
                                 final DenseMatrix64F gradient) {
        throw new UnsupportedOperationException("rootGradient not yet ported in local CDI prototype");
    }
}

final class OUActualizationStrategyGeneral extends AbstractOUActualizationStrategy {

    OUActualizationStrategyGeneral(final OUDelegateContext context) {
        super(context);
    }

    @Override
    public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                               final int precisionIndex,
                                               final double[] basisD,
                                               final double[] basisRotations) {
        integrator.setDiffusionStationaryVariance(precisionIndex, basisD, basisRotations);
    }

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                         final int precisionIndex,
                                                         final int[] probabilityIndices,
                                                         final double[] edgeLengths,
                                                         final double[] optimalRates,
                                                         final double[] basisD,
                                                         final double[] basisRotations,
                                                         final int updateCount) {
        integrator.updateOrnsteinUhlenbeckDiffusionMatrices(
                precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisD, basisRotations, updateCount);
    }

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                   final int precisionIndex,
                                                                   final int[] probabilityIndices,
                                                                   final double[] edgeLengths,
                                                                   final double[] optimalRates,
                                                                   final double[] basisD,
                                                                   final double[] basisRotations,
                                                                   final int updateCount) {
        integrator.updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
                precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisD, basisRotations, updateCount);
    }
}

final class OUActualizationStrategyDiagonal extends AbstractOUActualizationStrategy {

    OUActualizationStrategyDiagonal(final OUDelegateContext context) {
        super(context);
    }

    @Override
    public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                               final int precisionIndex,
                                               final double[] basisD,
                                               final double[] basisRotations) {
        integrator.setDiffusionStationaryVariance(precisionIndex, basisD, basisRotations);
    }

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                         final int precisionIndex,
                                                         final int[] probabilityIndices,
                                                         final double[] edgeLengths,
                                                         final double[] optimalRates,
                                                         final double[] basisD,
                                                         final double[] basisRotations,
                                                         final int updateCount) {
        integrator.updateOrnsteinUhlenbeckDiffusionMatrices(
                precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisD, basisRotations, updateCount);
    }

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                   final int precisionIndex,
                                                                   final int[] probabilityIndices,
                                                                   final double[] edgeLengths,
                                                                   final double[] optimalRates,
                                                                   final double[] basisD,
                                                                   final double[] basisRotations,
                                                                   final int updateCount) {
        throw new UnsupportedOperationException("Integrated diagonal OU strategy not yet ported in local CDI prototype");
    }
}

final class OUActualizationStrategyBlockDecomposition extends AbstractOUActualizationStrategy {

    OUActualizationStrategyBlockDecomposition(final OUDelegateContext context) {
        super(context);
    }

    @Override
    public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                               final int precisionIndex,
                                               final double[] basisD,
                                               final double[] basisRotations) {
        integrator.setDiffusionStationaryVariance(precisionIndex, basisD, basisRotations);
    }

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                         final int precisionIndex,
                                                         final int[] probabilityIndices,
                                                         final double[] edgeLengths,
                                                         final double[] optimalRates,
                                                         final double[] basisD,
                                                         final double[] basisRotations,
                                                         final int updateCount) {
        integrator.updateOrnsteinUhlenbeckDiffusionMatrices(
                precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisD, basisRotations, updateCount);
    }

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                   final int precisionIndex,
                                                                   final int[] probabilityIndices,
                                                                   final double[] edgeLengths,
                                                                   final double[] optimalRates,
                                                                   final double[] basisD,
                                                                   final double[] basisRotations,
                                                                   final int updateCount) {
        throw new UnsupportedOperationException("Integrated block OU strategy not yet ported in local CDI prototype");
    }
}
