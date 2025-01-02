package dr.math.geodesics;

public interface Manifold {
    void projectTangent(double[] tangent, double[] point);

    void projectPoint(double[] point);
    void geodesic(double[] point, double[] velocity, double t);
}
