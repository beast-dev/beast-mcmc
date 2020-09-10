package dr.inference.distribution;

public class NormalStatisticsHelper {

    public interface NormalStatisticsProvider {

        double getNormalMean();

        double getNormalPrecision();
    }

    public interface IndependentNormalStatisticsProvider {

        double getNormalMean(int dim);

        double getNormalPrecision(int dim);

    }

    public interface NormalMatrixStatisticsProvider {

        double getScalarPrecision(); //TODO: this should probably be removed

        double[] precisionMeanProduct(int col);

        double[][] getColPrecision(int col);

    }

}
