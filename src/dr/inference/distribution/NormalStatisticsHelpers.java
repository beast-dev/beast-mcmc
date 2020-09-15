package dr.inference.distribution;

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

public class NormalStatisticsHelpers {

    public interface NormalStatisticsHelper {

        MatrixNormalStatisticsHelper matrixNormalHelper(int nRows, int nCols);

    }

    public interface MatrixNormalStatisticsProvider extends NormalStatisticsHelper {

        double getNormalMean(int row, int col);

        double[] getColumnMean(int col);

        double[][] getColumnPrecision(int col);

        @Override
        default MatrixNormalStatisticsHelper matrixNormalHelper(int nRows, int nCols) {
            return new MatrixToMatrixAdaptor(this);
        }

    }

    public interface IndependentNormalStatisticsProvider extends NormalStatisticsHelper {

        double getNormalMean(int dim);

        double getNormalPrecision(int dim);

        @Override
        default MatrixNormalStatisticsHelper matrixNormalHelper(int nRows, int nCols) {
            return new IndependentToMatrixAdaptor(this, nRows, nCols);
        }

    }

    public interface MatrixNormalStatisticsHelper {

        double getScalarPrecision(); //TODO: this should probably be removed

        double[] precisionMeanProduct(int col);

        double[][] getColumnPrecision(int col); // TODO: remove

    }

    public static class MatrixToMatrixAdaptor implements MatrixNormalStatisticsHelper {
        private final MatrixNormalStatisticsProvider provider;

        public MatrixToMatrixAdaptor(MatrixNormalStatisticsProvider provider) {
            this.provider = provider;
        }


        @Override
        public double getScalarPrecision() {
            return provider.getColumnPrecision(0)[0][0]; //TODO: this probably isn't right
        }

        @Override
        public double[] precisionMeanProduct(int col) {

            double[] mean = provider.getColumnMean(col);
            double[][] prec = provider.getColumnPrecision(col);

            Vector meanVec = new Vector(mean);
            Matrix precMat = new Matrix(prec);

            try {
                Vector product = precMat.product(meanVec);
                return product.toComponents();
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            return null;
        }

        @Override
        public double[][] getColumnPrecision(int col) {
            return provider.getColumnPrecision(col);
        }
    }


    public static class IndependentToMatrixAdaptor implements MatrixNormalStatisticsProvider, MatrixNormalStatisticsHelper {
        private final IndependentNormalStatisticsProvider statistcisProvider;
        private final int nRows;
        private final int nCols;

        public IndependentToMatrixAdaptor(IndependentNormalStatisticsProvider statistcisProvider, int nRows, int nCols) {
            this.statistcisProvider = statistcisProvider;
            this.nRows = nRows;
            this.nCols = nCols;
        }

        @Override
        public double getNormalMean(int row, int col) {
            return statistcisProvider.getNormalMean(col * nRows + row); //TODO: check that indexing is correct
        }

        @Override
        public double[] getColumnMean(int col) {
            double[] mean = new double[nRows];

            for (int i = 0; i < nRows; i++) {
                mean[i] = getNormalMean(i, col);
            }

            return mean;
        }

        public double getColumnPrecisionDiagonal(int row, int col) {
            return statistcisProvider.getNormalPrecision(col + nCols * row);
        }

        @Override
        public double[][] getColumnPrecision(int col) {
            double[][] prec = new double[nRows][nRows];

            for (int i = 0; i < nRows; i++) {
                double p = getColumnPrecisionDiagonal(i, col); //TODO: check that indexing is correct
                prec[i][i] = p;
            }

            return prec;
        }

        @Override
        public double getScalarPrecision() {
            return statistcisProvider.getNormalPrecision(0);
        }

        @Override
        public double[] precisionMeanProduct(int col) {
            double[] prod = new double[nRows];
            for (int row = 0; row < nRows; row++) {
                prod[row] = getColumnPrecisionDiagonal(row, col) * getNormalMean(row, col);
            }
            return prod;
        }

        @Override
        public MatrixNormalStatisticsHelper matrixNormalHelper(int nRows, int nCols) {
            return this;
        }

    }

}
