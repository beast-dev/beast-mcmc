package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

public final class GaussianCanonicalTipObservation implements CanonicalTipObservationModel {

    private final int latentDimension;
    private final int observationDimension;
    private final double[] observation;
    private final double[] linkMatrix;
    private final double[] offset;
    private final double[] covariance;

    public GaussianCanonicalTipObservation(final int latentDimension,
                                           final int observationDimension,
                                           final double[] observation,
                                           final double[] linkMatrix,
                                           final double[] offset,
                                           final double[] covariance) {
        if (latentDimension < 1 || observationDimension < 1) {
            throw new IllegalArgumentException("dimensions must be positive");
        }
        requireLength(observation, observationDimension, "observation");
        requireLength(linkMatrix, observationDimension * latentDimension, "linkMatrix");
        requireLength(offset, observationDimension, "offset");
        requireLength(covariance, observationDimension * observationDimension, "covariance");
        this.latentDimension = latentDimension;
        this.observationDimension = observationDimension;
        this.observation = observation.clone();
        this.linkMatrix = linkMatrix.clone();
        this.offset = offset.clone();
        this.covariance = covariance.clone();
    }

    @Override
    public int getLatentDimension() {
        return latentDimension;
    }

    @Override
    public int getObservationDimension() {
        return observationDimension;
    }

    @Override
    public TipObservationMode getMode() {
        return TipObservationMode.GAUSSIAN_LINK;
    }

    @Override
    public boolean isEmpty() {
        return false;
    }

    @Override
    public void fillChildCanonicalState(final CanonicalGaussianState out,
                                        final TipObservationModelWorkspace workspace) {
        if (out.getDimension() != latentDimension) {
            throw new IllegalArgumentException("Output state dimension mismatch");
        }
        workspace.ensureCapacity(latentDimension, observationDimension);

        System.arraycopy(covariance, 0, workspace.covarianceScratch, 0,
                observationDimension * observationDimension);
        MatrixOps.invertSPD(workspace.covarianceScratch, workspace.precision, observationDimension);

        for (int i = 0; i < observationDimension; ++i) {
            workspace.shiftedObservation[i] = observation[i] - offset[i];
        }
        MatrixOps.matVec(
                workspace.precision,
                workspace.shiftedObservation,
                workspace.precisionTimesShifted,
                observationDimension);

        transposeObservationByLatent(linkMatrix, workspace.linkTranspose);
        MatrixOps.matMul(
                workspace.linkTranspose,
                workspace.precision,
                workspace.tempLatentByObservation,
                latentDimension,
                observationDimension,
                observationDimension);
        MatrixOps.matMul(
                workspace.tempLatentByObservation,
                linkMatrix,
                out.precision,
                latentDimension,
                observationDimension,
                latentDimension);
        MatrixOps.matVec(
                workspace.linkTranspose,
                workspace.precisionTimesShifted,
                out.information,
                latentDimension,
                observationDimension);

        final double logDet = MatrixOps.logDeterminant(covariance, observationDimension);
        double quadratic = 0.0;
        for (int i = 0; i < observationDimension; ++i) {
            quadratic += workspace.shiftedObservation[i] * workspace.precisionTimesShifted[i];
        }
        out.logNormalizer = 0.5 * (observationDimension * MatrixOps.LOG_TWO_PI + logDet + quadratic);
    }

    @Override
    public void fillParentMessage(final CanonicalGaussianTransition transition,
                                  final CanonicalTransitionMomentProvider momentProvider,
                                  final double branchLength,
                                  final TipParentMessageWorkspace workspace,
                                  final CanonicalGaussianMessageOps.Workspace gaussianWorkspace,
                                  final CanonicalGaussianState out) {
        fillChildCanonicalState(workspace.childStateScratch, workspace.observationModelWorkspace);
        CanonicalGaussianMessageOps.pushBackward(
                workspace.childStateScratch, transition, gaussianWorkspace, out);
    }

    @Override
    public CanonicalTipObservationModel copy() {
        return new GaussianCanonicalTipObservation(
                latentDimension,
                observationDimension,
                observation,
                linkMatrix,
                offset,
                covariance);
    }

    private void transposeObservationByLatent(final double[] source, final double[] target) {
        for (int obs = 0; obs < observationDimension; ++obs) {
            final int obsOffset = obs * latentDimension;
            for (int latent = 0; latent < latentDimension; ++latent) {
                target[latent * observationDimension + obs] = source[obsOffset + latent];
            }
        }
    }

    private static void requireLength(final double[] array, final int expected, final String name) {
        if (array == null || array.length != expected) {
            throw new IllegalArgumentException(name + " length must be " + expected);
        }
    }
}
