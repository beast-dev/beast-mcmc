package dr.evomodel.treedatalikelihood.continuous.observationmodel;

public final class TipObservationModelWorkspace {

    double[] shiftedObservation = new double[0];
    double[] precision = new double[0];
    double[] precisionTimesShifted = new double[0];
    double[] linkTranspose = new double[0];
    double[] tempLatentByObservation = new double[0];
    double[] covarianceScratch = new double[0];

    public void ensureCapacity(final int latentDimension, final int observationDimension) {
        if (shiftedObservation.length < observationDimension) {
            shiftedObservation = new double[observationDimension];
            precisionTimesShifted = new double[observationDimension];
        }
        final int observationSquare = observationDimension * observationDimension;
        if (precision.length < observationSquare) {
            precision = new double[observationSquare];
            covarianceScratch = new double[observationSquare];
        }
        final int latentByObservation = latentDimension * observationDimension;
        if (linkTranspose.length < latentByObservation) {
            linkTranspose = new double[latentByObservation];
            tempLatentByObservation = new double[latentByObservation];
        }
    }
}
