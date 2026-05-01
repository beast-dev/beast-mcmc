package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

public final class PartialIdentityTipObservationWorkspace {

    public final TipObservationPartition partition;
    public final double[] reducedPrecision;
    public final double[] reducedCovariance;
    public final double[] reducedInformation;
    public final double[] reducedMean;
    public final double[] missingMean;
    public final MatrixOps.MissingInversionWorkspace inversionWorkspace;

    public PartialIdentityTipObservationWorkspace(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.partition = new TipObservationPartition(dimension);
        this.reducedPrecision = new double[4 * dimension * dimension];
        this.reducedCovariance = new double[4 * dimension * dimension];
        this.reducedInformation = new double[2 * dimension];
        this.reducedMean = new double[2 * dimension];
        this.missingMean = new double[dimension];
        this.inversionWorkspace = new MatrixOps.MissingInversionWorkspace(2 * dimension);
    }
}
