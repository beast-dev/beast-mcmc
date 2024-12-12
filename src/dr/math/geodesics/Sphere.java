package dr.math.geodesics;

public class Sphere implements Manifold {


    public Sphere() {

    }

    @Override
    public void projectTangent(double[] tangent, double[] point) {
        int dim = point.length;

        double radiusSquared = 0;
        for (int i = 0; i < dim; i++) {
            radiusSquared += point[i] * point[i];
        }


        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                tangent[i] -= point[j] * point[i] * tangent[j] / radiusSquared;
            }
        }
    }

    @Override
    public void geodesic(double[] point, double[] velocity, double t) {
        // assumes velocity is already orthogonal to point

        int dim = point.length;
        double alpha = 0;
        for (int i = 0; i < dim; i++) {
            alpha += velocity[i] * velocity[i];
        }
        alpha = Math.sqrt(alpha);

        double cat = Math.cos(alpha * t);
        double sat = Math.sin(alpha * t);

        for (int i = 0; i < dim; i++) {
            point[i] = point[i] * cat + velocity[i] * sat / alpha;
            velocity[i] = velocity[i] * cat - point[i] * sat * alpha;
        }
    }


}
