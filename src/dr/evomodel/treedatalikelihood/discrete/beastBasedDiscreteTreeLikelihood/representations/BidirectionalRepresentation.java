package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

public interface BidirectionalRepresentation
        extends PostOrderRepresentation, PreOrderRepresentation {

    @Override
    default boolean storesPartialsInStandardBasis() {
        return PostOrderRepresentation.super.storesPartialsInStandardBasis();
    }
}
