package dr.inference.model;

import dr.app.bss.Utils;
import dr.evomodel.substmodel.ColtEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.compoundSymmetricMatrix;

public interface BoundedSpace {

    boolean isWithinBounds(double[] values);

    IntersectionDistances distancesToBoundary(double[] origin, double[] direction);

    default double forwardDistanceToBoundary(double[] origin, double[] direction) {
        return distancesToBoundary(origin, direction).forwardDistance;
    }

    class IntersectionDistances {
        public final double forwardDistance;
        public final double backwardDistance;
        public final double totalDistance;

        public IntersectionDistances(double forwardDistance, double backwardDistance) {
            this.forwardDistance = forwardDistance;
            this.backwardDistance = backwardDistance;
            this.totalDistance = forwardDistance + backwardDistance;
        }
    }


    class Correlation implements BoundedSpace {

        private static final boolean DEBUG = false;
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

        @Override
        public IntersectionDistances distancesToBoundary(double[] origin, double[] direction) {
            double[] x = new double[origin.length];
            System.arraycopy(direction, 0, x, 0, x.length);

            SymmetricMatrix Y = compoundCorrelationSymmetricMatrix(origin, dim);
            SymmetricMatrix X = compoundSymmetricMatrix(0.0, x, dim);

//        SymmetricMatrix Xinv = X.inverse();
            SymmetricMatrix Yinv = Y.inverse();
            final Matrix Z;

            try {
                Z = Yinv.product(X);
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("illegal dimensions");
            }

            ColtEigenSystem eigenSystem = new ColtEigenSystem(dim);
            EigenDecomposition decomposition = eigenSystem.decomposeMatrix(Z.toComponents()); //TODO: only need smallest magnitude eigenvalues
            double[] values = decomposition.getEigenValues();

            double maxNegative = 0;
            double maxPositive = 0;
            for (int i = 0; i < values.length; i++) {
                double value = values[i];
                if (value < 0 && value < maxNegative) {
                    maxNegative = value;
                } else if (value >= 0 & value > maxPositive) {
                    maxPositive = value;
                }
            }

            if (DEBUG) {
                System.out.print("Eigenvalues: ");
                Utils.printArray(values);

                Matrix S = new SymmetricMatrix(dim, dim);
                Matrix T = new SymmetricMatrix(dim, dim);
                for (int i = 0; i < dim; i++) {
                    S.set(i, i, 1);
                    T.set(i, i, 1);
                    for (int j = (i + 1); j < dim; j++) {
                        double y = Y.toComponents()[i][j];
                        double z = X.toComponents()[i][j];
                        double valueS = y - z / maxNegative;
                        double valueT = y - z / maxPositive;
                        S.set(i, j, valueS);
                        S.set(j, i, valueS);
                        T.set(i, j, valueT);
                        T.set(j, i, valueT);
                    }
                }
                try {
                    System.out.println("neg: \n\tt = " + maxNegative);
                    System.out.println("\tdet = " + S.determinant());
                    System.out.println(S);
                    System.out.println("pos: \n\tt = " + maxPositive);
                    System.out.println("\tdet = " + T.determinant());
                    System.out.println(T);
                } catch (IllegalDimension illegalDimension) {
                    illegalDimension.printStackTrace();
                }
            }

            return new IntersectionDistances(1 / maxPositive, -1 / maxNegative);
        }

    }

}
