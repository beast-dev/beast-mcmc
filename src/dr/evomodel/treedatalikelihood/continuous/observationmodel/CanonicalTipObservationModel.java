package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

public interface CanonicalTipObservationModel {

    int getLatentDimension();

    int getObservationDimension();

    TipObservationMode getMode();

    boolean isEmpty();

    void fillChildCanonicalState(CanonicalGaussianState out,
                                 TipObservationModelWorkspace workspace);

    CanonicalTipObservationModel copy();
}
