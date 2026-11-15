package dr.evomodel.treedatalikelihood.continuousNew.observationmodel;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;

import java.util.Arrays;

/**
 * Exact linear observation of a latent child state, {@code z = H y}.
 *
 * <p>This is used by the homogeneous block-basis OU path for partially
 * observed natural-coordinate tips, where {@code H} is dense in the latent
 * block basis. It is not a finite-noise Gaussian link and therefore cannot be
 * represented as a child-side canonical state.</p>
 */
public final class LinearExactTipObservationModel implements CanonicalTipObservationModel {

    private final int latentDimension;
    private final int observationDimension;
    private final int[] observedIndices;
    private final double[] shiftedObservation;
    private final double[] linkMatrix;

    public LinearExactTipObservationModel(final int latentDimension,
                                          final int observationDimension,
                                          final int[] observedIndices,
                                          final double[] shiftedObservation,
                                          final double[] linkMatrix) {
        if (latentDimension < 1) {
            throw new IllegalArgumentException("latentDimension must be positive");
        }
        if (observationDimension < 1 || observationDimension >= latentDimension) {
            throw new IllegalArgumentException(
                    "Linear exact observations are for genuine partial tips: "
                            + observationDimension + " vs " + latentDimension);
        }
        requireLength(observedIndices, observationDimension, "observedIndices");
        requireLength(shiftedObservation, observationDimension, "shiftedObservation");
        requireLength(linkMatrix, observationDimension * latentDimension, "linkMatrix");
        this.latentDimension = latentDimension;
        this.observationDimension = observationDimension;
        this.observedIndices = Arrays.copyOf(observedIndices, observationDimension);
        this.shiftedObservation = Arrays.copyOf(shiftedObservation, observationDimension);
        this.linkMatrix = Arrays.copyOf(linkMatrix, observationDimension * latentDimension);
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
        return TipObservationMode.PARTIAL_LINEAR_EXACT;
    }

    @Override
    public boolean isEmpty() {
        return false;
    }

    public int observedIndexAt(final int observationIndex) {
        return observedIndices[observationIndex];
    }

    public double shiftedValueAt(final int observationIndex) {
        return shiftedObservation[observationIndex];
    }

    public double linkValueAt(final int observationIndex, final int latentIndex) {
        return linkMatrix[observationIndex * latentDimension + latentIndex];
    }

    public void copyFrom(final LinearExactTipObservationModel source) {
        if (source.latentDimension != latentDimension
                || source.observationDimension != observationDimension) {
            throw new IllegalArgumentException(
                    "Linear exact observation shape mismatch: "
                            + source.observationDimension + "x" + source.latentDimension
                            + " vs " + observationDimension + "x" + latentDimension);
        }
        System.arraycopy(source.observedIndices, 0, observedIndices, 0, observationDimension);
        System.arraycopy(source.shiftedObservation, 0, shiftedObservation, 0, observationDimension);
        System.arraycopy(source.linkMatrix, 0, linkMatrix, 0,
                observationDimension * latentDimension);
    }

    @Override
    public void fillParentMessage(final CanonicalGaussianTransition transition,
                                  final int childNodeIndex,
                                  final CanonicalTransitionMomentProvider momentProvider,
                                  final double branchLength,
                                  final TipParentMessageWorkspace workspace,
                                  final CanonicalGaussianMessageOps.Workspace gaussianWorkspace,
                                  final CanonicalGaussianState out) {
        if (momentProvider == null) {
            throw new UnsupportedOperationException(
                    "Linear exact observations require a CanonicalTransitionMomentProvider.");
        }
        workspace.linearExactProjection.projectObservedLinearChildToParent(
                this, childNodeIndex, momentProvider, branchLength, out);
    }

    @Override
    public void fillChildCanonicalState(final CanonicalGaussianState out,
                                        final TipObservationModelWorkspace workspace) {
        throw new UnsupportedOperationException(
                "Exact linear observations cannot be represented as finite canonical child states; "
                        + "use fillParentMessage instead.");
    }

    @Override
    public CanonicalTipObservationModel copy() {
        return new LinearExactTipObservationModel(
                latentDimension,
                observationDimension,
                observedIndices,
                shiftedObservation,
                linkMatrix);
    }

    private static void requireLength(final int[] array, final int expected, final String name) {
        if (array == null || array.length < expected) {
            throw new IllegalArgumentException(name + " length must be at least " + expected);
        }
    }

    private static void requireLength(final double[] array, final int expected, final String name) {
        if (array == null || array.length < expected) {
            throw new IllegalArgumentException(name + " length must be at least " + expected);
        }
    }
}
