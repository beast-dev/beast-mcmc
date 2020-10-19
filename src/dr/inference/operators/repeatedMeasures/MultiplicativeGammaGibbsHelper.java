package dr.inference.operators.repeatedMeasures;

public interface MultiplicativeGammaGibbsHelper {

    double computeSumSquaredErrors(int column);

    int getRowDimension();

    int getColumnDimension();

    public class NormalizedMultiplicativeGammaGibbsHelper implements MultiplicativeGammaGibbsHelper {
        private final int rowDimension;
        private final int colDimension;

        public NormalizedMultiplicativeGammaGibbsHelper(int rowDimension, int colDimension) {
            this.rowDimension = rowDimension;
            this.colDimension = colDimension;

        }

        @Override
        public double computeSumSquaredErrors(int column) {
            return 1.0;
        }

        @Override
        public int getRowDimension() {
            return rowDimension;
        }

        @Override
        public int getColumnDimension() {
            return colDimension;
        }
    }
}
