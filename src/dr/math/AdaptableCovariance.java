package dr.math;

import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Xiang Ji
 */
public class AdaptableCovariance {

    final private int dim;
    final private double[][] empirical;
    final private AdaptableVector.Default means;

    protected int updates;

    public AdaptableCovariance(int dim) {
        this.dim = dim;
        this.empirical = new double[dim][dim];
        this.means = new AdaptableVector.Default(dim);

        updates = 0;
        for (int i = 0; i < dim; i++) {
            empirical[i][i] = 1.0;
        }
    }

    public int getUpdateCount() { return updates; }

    public void update(ReadableVector x) {

        assert (x.getDim() == dim);

        if (shouldUpdate()) {

            ++updates;

            means.update(x);

            if (updates > 1) {
                updateVariance(x);
            }
        }
    }

    public ReadableMatrix getCovariance() {
        double[][] copyEmpirical = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            copyEmpirical[i] = empirical[i].clone();
        }
        return new WrappedMatrix.ArrayOfArray(copyEmpirical);
    }

    public ReadableVector getMean() {
        return means.getMean();
    }

    protected boolean shouldUpdate() { return true; }

    private void updateVariance(ReadableVector x) {
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                empirical[i][j] = calculateCovariance(empirical[i][j], x, i, j);
                empirical[j][i] = empirical[i][j];
            }
        }
    }

    private double calculateCovariance(double currentMatrixEntry, ReadableVector values, int firstIndex, int secondIndex) {

        double result = currentMatrixEntry * (updates - 2);
        result += (values.get(firstIndex) * values.get(secondIndex));
        result += ((updates - 1) * means.getOldMeans(firstIndex) * means.getOldMeans(secondIndex) - updates * means.getNewMeans(firstIndex) * means.getNewMeans(secondIndex));
        result /= ((double)(updates - 1));

        return result;
    }

    public static class WithSubsampling extends AdaptableCovariance {

        final int maxUpdates = 1000;

        public WithSubsampling(int dim) {
            super(dim);
        }

        @Override
        protected boolean shouldUpdate() {
             // TODO Add logic in subclass to control how often updates are made
//             return true;
            return updates < maxUpdates;
         }

    }
}
