package dr.math.distributions;

public interface ConvexSpaceRandomGenerator extends RandomGenerator {

    LineThroughPoints distanceToEdge(double[] origin, double[] draw);

    boolean isUniform();

    class LineThroughPoints {
        public final double forwardDistance;
        public final double backwardDistance;
        public final double totalDistance;

        public LineThroughPoints(double forwardDistance, double backwardDistance) {
            this.forwardDistance = forwardDistance;
            this.backwardDistance = backwardDistance;
            this.totalDistance = forwardDistance + backwardDistance;
        }
    }

}
