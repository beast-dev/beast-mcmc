package dr.inference.model;

import org.ejml.alg.dense.decomposition.chol.CholeskyDecompositionCommon_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

public interface BoundedSpace {

    boolean isWithinBounds(double[] values);


    class Correlation implements BoundedSpace {

        private final int dim;

        public Correlation(int dim) {
            this.dim = dim;
        }


        @Override
        public boolean isWithinBounds(double[] x) {

            DenseMatrix64F C;
            double[] values = new double[x.length];
            System.arraycopy(x, 0, values, 0, x.length);

            if (values.length == dim * dim) {
                C = DenseMatrix64F.wrap(dim, dim, values);
                for (int i = 0; i < dim; i++) {
                    if (C.get(i, i) != 1.0) {
                        return false;
                    }
                }

            } else if (values.length == dim * (dim - 1) / 2) {
                int ind = 0;
                C = new DenseMatrix64F(dim, dim);
                for (int i = 0; i < dim; i++) {
                    C.set(i, i, 1.0);
                    for (int j = i + 1; j < dim; j++) {
                        C.set(i, j, values[ind]);
                        C.set(j, i, values[ind]);
                        ind++;
                    }
                }
            } else {
                throw new RuntimeException("incompatible dimensions");
            }


            CholeskyDecomposition<DenseMatrix64F> chol = DecompositionFactory.chol(dim, true);
            boolean isDecomposable = chol.decompose(C); // in place decomposition
            if (!isDecomposable) {
                return false;
            }

            for (int i = 0; i < dim; i++) {
                if (C.get(i, i) <= 0) {
                    return false;
                }
            }

            return true;
        }

    }

}
