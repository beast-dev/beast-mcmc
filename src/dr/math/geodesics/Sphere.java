package dr.math.geodesics;

import dr.math.MathUtils;

public class Sphere implements Manifold {

    private double radius;

    public Sphere(double radius) {
        this.radius = radius;
    }

    public Sphere() {
        this(1);
    }

    @Override
    public void projectTangent(double[] tangent, double[] point) {
        int dim = point.length;
        double[] originalTangent = new double[dim];
        System.arraycopy(tangent, 0, originalTangent, 0, dim);

        double radiusSquared = radius * radius;


        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                tangent[i] -= point[j] * point[i] * originalTangent[j] / radiusSquared;
            }
        }
    }

    @Override
    public void projectPoint(double[] point) {
        double currentRadius = MathUtils.getL2Norm(point);
        double ratio = radius / currentRadius;

        for (int i = 0; i < point.length; i++) {
            point[i] = point[i] * ratio;
        }
    }

    @Override
    public void geodesic(double[] point, double[] velocity, double t) {
        // assumes velocity is already orthogonal to point

        int dim = point.length;
        double alpha = MathUtils.getL2Norm(velocity) / radius;

        double cat = Math.cos(alpha * t);
        double sat = Math.sin(alpha * t);

        for (int i = 0; i < dim; i++) {
            point[i] = point[i] * cat + velocity[i] * sat / alpha;
            velocity[i] = velocity[i] * cat - point[i] * sat * alpha;
        }
    }

    @Override
    public void initialize(double[] values) {
        radius = MathUtils.getL2Norm(values);
    }


}
