package dr.inference.timeseries.representation;

/**
 * Bridge from semantic process definitions to computational representations.
 */
public interface RepresentableProcess {

    <T> boolean supportsRepresentation(Class<T> representationClass);

    <T> T getRepresentation(Class<T> representationClass);
}
