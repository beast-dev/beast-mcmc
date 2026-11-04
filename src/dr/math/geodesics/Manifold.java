package dr.math.geodesics;

import dr.math.matrixAlgebra.WrappedMatrix;

public interface Manifold {
    void projectTangent(double[] tangent, double[] point);

    void projectPoint(double[] point);

    void geodesic(double[] point, double[] velocity, double t);

    void initialize(double[] values);

    default void initialize(WrappedMatrix matrix) {
        throw new RuntimeException("Only valid for matrix manifolds");
    }

    default boolean isMatrix() {
        return false;
    }

    public interface MatrixManifold extends Manifold {
        default boolean isMatrix() {
            return true;
        }

        default void initialize(double[] values) {
            throw new RuntimeException("Not valid for matrix manifolds");
        }

        void initialize(WrappedMatrix matrix);
    }

}
