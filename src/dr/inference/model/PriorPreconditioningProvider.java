package dr.inference.model;

public interface PriorPreconditioningProvider {
    double getStandardDeviation(int index);
    int getDimension();
}
