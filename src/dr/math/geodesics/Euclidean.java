package dr.math.geodesics;

public class Euclidean implements Manifold {
    @Override
    public void projectTangent(double[] tangent, double[] point) {
        // do nothing
    }

    @Override
    public void projectPoint(double[] point) {
        // do nothing
    }

    @Override
    public void geodesic(double[] point, double[] velocity, double t) {
        for (int i = 0; i < point.length; i++) {
            point[i] = point[i] + t * velocity[i];
        }
    }

    @Override
    public void initialize(double[] values) {
        // do nothing
    }
}
