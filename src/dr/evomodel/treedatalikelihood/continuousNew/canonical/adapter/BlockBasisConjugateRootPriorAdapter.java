package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evomodel.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;

/**
 * Conjugate root prior expressed in the homogeneous block basis
 * {@code z = R^{-1}(y - beta)}.
 */
final class BlockBasisConjugateRootPriorAdapter implements CanonicalRootPrior {

    private final ConjugateRootTraitPrior prior;
    private final Parameter stationaryMean;
    private final BlockBasisGlobalSnapshot snapshot;
    private final int dim;

    private final double[] meanZ;
    private final double[] posteriorMeanZ;
    private final double[] scratchCovarianceFlat;
    private final double[] rootGradientZ;
    private final CanonicalGaussianState scratchPriorState;
    private final CanonicalGaussianState scratchCombinedState;
    private final CanonicalGaussianMessageOps.Workspace workspace;
    private final GaussianFormConverter.Workspace converterWorkspace;

    BlockBasisConjugateRootPriorAdapter(final ConjugateRootTraitPrior prior,
                                        final Parameter stationaryMean,
                                        final BlockBasisGlobalSnapshot snapshot) {
        this.prior = prior;
        this.stationaryMean = stationaryMean;
        this.snapshot = snapshot;
        this.dim = snapshot.getDimension();
        this.meanZ = new double[dim];
        this.posteriorMeanZ = new double[dim];
        this.scratchCovarianceFlat = new double[dim * dim];
        this.rootGradientZ = new double[dim];
        this.scratchPriorState = new CanonicalGaussianState(dim);
        this.scratchCombinedState = new CanonicalGaussianState(dim);
        this.workspace = new CanonicalGaussianMessageOps.Workspace(dim);
        this.converterWorkspace = new GaussianFormConverter.Workspace();
        this.converterWorkspace.ensureDim(dim);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public void fillRootPriorState(final double[] traitCovariance,
                                   final CanonicalGaussianState out) {
        final double kappa0 = prior.getPseudoObservations();
        if (!(kappa0 > 0.0) || Double.isInfinite(kappa0)) {
            throw new UnsupportedOperationException(
                    "Block-basis canonical root seeding requires finite positive pseudo-observations.");
        }

        fillRootMeanZ();
        final double invKappa0 = 1.0 / kappa0;
        for (int k = 0; k < dim * dim; ++k) {
            scratchCovarianceFlat[k] = traitCovariance[k] * invKappa0;
        }
        GaussianFormConverter.fillStateFromMoments(
                meanZ,
                scratchCovarianceFlat,
                dim,
                converterWorkspace,
                out);
    }

    @Override
    public double computeLogMarginalLikelihood(final CanonicalGaussianState rootMessage,
                                               final double[] traitCovariance) {
        if (isFixedRoot()) {
            fillRootMeanZ();
            double quadratic = 0.0;
            double linear = 0.0;
            for (int i = 0; i < dim; ++i) {
                linear += rootMessage.information[i] * meanZ[i];
                for (int j = 0; j < dim; ++j) {
                    quadratic += meanZ[i] * rootMessage.precision[i * dim + j] * meanZ[j];
                }
            }
            return -0.5 * quadratic + linear - rootMessage.logNormalizer;
        }

        fillRootPriorState(traitCovariance, scratchPriorState);
        CanonicalGaussianMessageOps.combineStates(rootMessage, scratchPriorState, scratchCombinedState);
        return CanonicalGaussianMessageOps.normalizationShift(scratchCombinedState, workspace);
    }

    @Override
    public boolean isFixedRoot() {
        return Double.isInfinite(prior.getPseudoObservations());
    }

    @Override
    public void fillFixedRootValue(final double[] out) {
        fillRootMeanZ();
        System.arraycopy(meanZ, 0, out, 0, dim);
    }

    @Override
    public void accumulateRootMeanGradient(final CanonicalGaussianState rootMessage,
                                           final double[] traitCovariance,
                                           final double[] gradOut) {
        if (isFixedRoot()) {
            fillRootMeanZ();
            for (int i = 0; i < dim; ++i) {
                double gradient = rootMessage.information[i];
                for (int j = 0; j < dim; ++j) {
                    gradient -= rootMessage.precision[i * dim + j] * meanZ[j];
                }
                rootGradientZ[i] = gradient;
            }
            accumulateNaturalRootMeanGradient(rootGradientZ, gradOut);
            return;
        }

        fillRootPriorState(traitCovariance, scratchPriorState);
        CanonicalGaussianMessageOps.combineStates(rootMessage, scratchPriorState, scratchCombinedState);
        GaussianFormConverter.fillMomentsFromState(
                scratchCombinedState,
                posteriorMeanZ,
                scratchCovarianceFlat,
                dim,
                converterWorkspace);

        fillRootMeanZ();
        for (int i = 0; i < dim; ++i) {
            double gradient = 0.0;
            for (int j = 0; j < dim; ++j) {
                gradient += scratchPriorState.precision[i * dim + j] * (posteriorMeanZ[j] - meanZ[j]);
            }
            rootGradientZ[i] = gradient;
        }
        accumulateNaturalRootMeanGradient(rootGradientZ, gradOut);
    }

    @Override
    public double getDiffusionScale() {
        final double kappa0 = prior.getPseudoObservations();
        if (!(kappa0 > 0.0) || Double.isInfinite(kappa0)) {
            return 0.0;
        }
        return 1.0 / kappa0;
    }

    private void fillRootMeanZ() {
        snapshot.ensureCurrent();
        final DenseMatrix64F t = snapshot.getRotationInverse();
        final double[] rootMean = prior.getMean();
        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dim; ++j) {
                sum += t.get(i, j) * (rootMean[j] - stationaryMean.getParameterValue(j));
            }
            meanZ[i] = sum;
        }
    }

    private void accumulateNaturalRootMeanGradient(final double[] gradientZ,
                                                   final double[] gradOut) {
        final DenseMatrix64F t = snapshot.getRotationInverse();
        for (int j = 0; j < dim; ++j) {
            double sum = 0.0;
            for (int i = 0; i < dim; ++i) {
                sum += t.get(i, j) * gradientZ[i];
            }
            gradOut[j] += sum;
        }
    }
}
