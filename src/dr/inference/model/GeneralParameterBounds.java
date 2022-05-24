package dr.inference.model;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public interface GeneralParameterBounds {

    boolean satisfiesBounds(Parameter parameter);


    class CorrelationParameterBounds implements GeneralParameterBounds {

        private final int dim;

        public CorrelationParameterBounds(int dim) {
            this.dim = dim;
        }


        @Override
        public boolean satisfiesBounds(Parameter parameter) {

            DenseMatrix64F C;
            double[] c = parameter.getParameterValues();

            if (c.length == dim * dim) {
                C = DenseMatrix64F.wrap(dim, dim, parameter.getParameterValues());
                for (int i = 0; i < dim; i++) {
                    if (C.get(i, i) != 1.0) {
                        return false;
                    }
                }

            } else if (c.length == dim * (dim - 1) / 2) {
                int ind = 0;
                C = new DenseMatrix64F(dim, dim);
                for (int i = 0; i < dim; i++) {
                    C.set(i, i, 1.0);
                    for (int j = i + 1; j < dim; j++) {
                        C.set(i, j, c[ind]);
                        C.set(j, i, c[ind]);
                        ind++;
                    }
                }
            } else {
                throw new RuntimeException("incompatible dimensions");
            }


            double det = CommonOps.det(C);
            return det >= 0; // already checked if diagonals were 1
        }

    }

}
