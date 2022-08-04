package dr.inference.model;

import dr.app.bss.Utils;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;
import static dr.math.matrixAlgebra.SymmetricMatrix.compoundSymmetricMatrix;

public interface BoundedSpace extends GeneralBoundsProvider {

    boolean isWithinBounds(double[] values);

    IntersectionDistances distancesToBoundary(double[] origin, double[] direction, boolean isAtBoundary) throws HamiltonianMonteCarloOperator.NumericInstabilityException;

    double[] getNormalVectorAtBoundary(double[] position);

    default double forwardDistanceToBoundary(double[] origin, double[] direction, boolean isAtBoundary) throws HamiltonianMonteCarloOperator.NumericInstabilityException {
        return distancesToBoundary(origin, direction, isAtBoundary).forwardDistance;
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
        private static final boolean CHECK_WITHIN_BOUNDS = false;
        private static final double TOL = 0;
        private static final double BOUNDARY_TOL = 1e-9;
        private final int dim;

        private final DenseMatrix64F C;
        private final DenseMatrix64F V;
        private final DenseMatrix64F Cinv;
        private final DenseMatrix64F CinvV;

        public Correlation(int dim) {
            this.dim = dim;
            this.V = new DenseMatrix64F(dim, dim);
            this.C = new DenseMatrix64F(dim, dim);
            for (int i = 0; i < dim; i++) {
                C.set(i, i, 1);
            }
            this.Cinv = new DenseMatrix64F(dim, dim);
            this.CinvV = new DenseMatrix64F(dim, dim);
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

        private double[] robustTrajectoryEigenValues(double[] origin, double[] direction) {
            double t = MathUtils.nextDouble();
            double[] newOrigin = new double[origin.length];
            for (int i = 0; i < origin.length; i++) {
                newOrigin[i] = origin[i] + t * direction[i];
            }

            double[] shiftedValues = trajectoryEigenvalues(newOrigin, direction);
            for (int i = 0; i < shiftedValues.length; i++) {
                shiftedValues[i] -= t;
            }
            return shiftedValues;
        }

        private double[] trajectoryEigenvalues(double[] origin, double[] direction) {
            double[] x = new double[origin.length];
            System.arraycopy(direction, 0, x, 0, x.length); //TODO: is this necessary?

            int ind = 0;
            for (int i = 0; i < dim; i++) {
                for (int j = (i + 1); j < dim; j++) {
                    C.set(i, j, origin[ind]);
                    C.set(j, i, origin[ind]);
                    V.set(i, j, direction[ind]);
                    V.set(j, i, direction[ind]);
                    ind++;
                }
            }

            CommonOps.invert(C, Cinv);
            CommonOps.mult(Cinv, V, CinvV);

            org.ejml.interfaces.decomposition.EigenDecomposition<DenseMatrix64F> factory = new DecompositionFactory().eig(dim, false, false);

            if (!factory.decompose(CinvV)) throw new RuntimeException("Eigen decomposition failed.");

            double[] allValues = new double[dim];
            int nReal = 0;
            for (int i = 0; i < dim; i++) {
                Complex64F ev = factory.getEigenvalue(i);
                if (ev.isReal()) {
                    allValues[nReal] = ev.real;
                    nReal++;
                }
            }

            double[] values = new double[nReal];
            System.arraycopy(allValues, 0, values, 0, nReal);

            if (DEBUG) {
                System.out.println("Raw matrix to decompose: ");
                System.out.println(CinvV);
                System.out.print("Raw eigenvalues: ");
                Utils.printArray(values);
            }
            for (int i = 0; i < values.length; i++) {
                values[i] = 1 / values[i];
            }

            return values;

        }

        @Override
        public IntersectionDistances distancesToBoundary(double[] origin, double[] direction, boolean isAtBoundary) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

            if (CHECK_WITHIN_BOUNDS && !isWithinBounds(origin)) { // don't automatically check that it's inside the boundary
                if (isAtBoundary) {
                    SymmetricMatrix C = compoundCorrelationSymmetricMatrix(origin, dim);
                    double det;
                    try {
                        det = C.determinant();
                    } catch (IllegalDimension illegalDimension) {
                        illegalDimension.printStackTrace();
                        throw new RuntimeException();
                    }

                    if (Math.abs(det) > BOUNDARY_TOL) {
                        System.out.println(det);
                        throw new HamiltonianMonteCarloOperator.NumericInstabilityException();
                    }
                } else {
                    SymmetricMatrix C = compoundCorrelationSymmetricMatrix(origin, dim);
                    System.out.println(C);
                    try {
                        System.out.println(C.determinant());
                    } catch (IllegalDimension illegalDimension) {
                        illegalDimension.printStackTrace();
                        throw new RuntimeException();
                    }

                    throw new HamiltonianMonteCarloOperator.NumericInstabilityException();
                }
            }


            double values[] = robustTrajectoryEigenValues(origin, direction);

            double minNegative = Double.NEGATIVE_INFINITY;
            double minNegative2 = Double.NEGATIVE_INFINITY;
            double minPositive = Double.POSITIVE_INFINITY;
            double minPositive2 = Double.POSITIVE_INFINITY;

            for (int i = 0; i < values.length; i++) {
                double value = values[i];
                if (value < -TOL && value > minNegative) {
                    minNegative2 = minNegative;
                    minNegative = value;
                } else if (value < -TOL && value > minNegative2) {
                    minNegative2 = value;
                } else if (value >= TOL & value < minPositive) {
                    minPositive2 = minPositive;
                    minPositive = value;
                } else if (value >= TOL && value < minPositive2) {
                    minPositive2 = value;
                }
            }

            if (isAtBoundary) {

                if (Math.abs(minNegative) < minPositive) {
                    if (Math.abs(minNegative) < BOUNDARY_TOL) {
                        minNegative = minNegative2;
                    } else {
                        throw new RuntimeException("isAtBoundary = true but does not appear to be near boundary");
                    }
                } else {
                    if (minPositive < BOUNDARY_TOL) {
                        minPositive = minPositive2;
                    } else {
                        throw new RuntimeException("isAtBoundary = true but does not appear to be near boundary");
                    }
                }

            }

            minPositive = -minPositive;
            minNegative = -minNegative;

            if (DEBUG) {
                SymmetricMatrix Y = compoundCorrelationSymmetricMatrix(origin, dim);
                SymmetricMatrix X = compoundSymmetricMatrix(0.0, direction, dim);
                System.out.print("Eigenvalues: ");
                Utils.printArray(values);

                Matrix S = new SymmetricMatrix(dim, dim);
                Matrix T = new SymmetricMatrix(dim, dim);
                double absMax = 0.0;
                for (int i = 0; i < dim; i++) {
                    S.set(i, i, 1);
                    T.set(i, i, 1);
                    for (int j = (i + 1); j < dim; j++) {
                        double y = Y.toComponents()[i][j];
                        double z = X.toComponents()[i][j];
                        double valueS = y + z * minNegative;
                        double valueT = y + z * minPositive;
                        if (Math.abs(valueS) > absMax) absMax = Math.abs(valueS);
//                        if (Math.abs(valueT) > absMax) absMax = Math.abs(valueT);

                        S.set(i, j, valueS);
                        S.set(j, i, valueS);
                        T.set(i, j, valueT);
                        T.set(j, i, valueT);
                    }
                }
                double detY;
                try {
                    System.out.println("starting position: ");
                    detY = Y.determinant();
                    System.out.println("\tdet = " + detY);
                    System.out.println(Y);
                    System.out.println("direction:");
                    System.out.println(X);
                    System.out.println();
                    System.out.println("neg: \n\tt = " + minNegative);
                    System.out.println("\tdet = " + S.determinant());
                    System.out.println(S);
                    System.out.println("pos: \n\tt = " + minPositive);
                    System.out.println("\tdet = " + T.determinant());
                    System.out.println(T);
                } catch (IllegalDimension illegalDimension) {
                    illegalDimension.printStackTrace();
                    throw new RuntimeException();
                }

            }

            return new IntersectionDistances(minNegative, minPositive);
        }

        @Override
        public double[] getNormalVectorAtBoundary(double[] position) {
            double[] c = compoundCorrelationSymmetricMatrix(position, dim).toArrayComponents();
            DenseMatrix64F C = DenseMatrix64F.wrap(dim, dim, c);
            DenseMatrix64F A;
            try {
                A = EJMLUtils.computeRobustAdjugate(C);
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
                throw new RuntimeException();
            }

            double[] normalVector = new double[position.length];
            int ind = 0;
            for (int i = 0; i < dim; i++) {
                for (int j = (i + 1); j < dim; j++) {
                    normalVector[ind] = A.get(i, j);
                    ind++;
                }
            }

            return normalVector;
        }

    }

}
