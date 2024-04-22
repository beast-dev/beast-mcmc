package dr.inference.distribution;

public interface GammaStatisticsProvider {
    double getShape(int dim);

    double getRate(int dim);
}
