package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.ConvexSpaceRandomGenerator;
import jebl.math.Random;

public class ConvexSpaceRandomWalkOperator extends AbstractAdaptableOperator {
    private double window;
    private final ConvexSpaceRandomGenerator generator;
    private final Parameter parameter;
    private final Parameter updateIndex;

    public static final String CONVEX_RW = "convexSpaceRandomWalkOperator";
    public static final String WINDOW_SIZE = "relativeWindowSize";

    public ConvexSpaceRandomWalkOperator(Parameter parameter, ConvexSpaceRandomGenerator generator,
                                         Parameter updateIndex,
                                         double window, double weight) {
        setWeight(weight);

        this.updateIndex = updateIndex;
        this.parameter = parameter;
        this.generator = generator;
        this.window = window;
    }


    @Override
    public double doOperation() {
//        double[] sample = (double[]) generator.nextRandom();
        double[] values = parameter.getParameterValues();

        double[] sample = new double[values.length];
        double sum = 0;
        for (int i = 0; i < values.length; i++) {
            if (updateIndex == null || updateIndex.getParameterValue(i) == 1) {
                sample[i] = MathUtils.nextGaussian();
                sum += sample[i] * sample[i];
            }
        }
        double norm = Math.sqrt(sum);
        for (int i = 0; i < values.length; i++) {
            sample[i] = sample[i] / norm;
        }

        ConvexSpaceRandomGenerator.LineThroughPoints distances = generator.distanceToEdge(values, sample);
//        double u1 = Random.nextDouble() * distances.forwardDistance;
//        for (int i = 0; i < values.length; i++) {
//            sample[i] = values[i] + (sample[i] - values[i]) * u1;
//        }


        double t = window * Random.nextDouble();
        t = RandomWalkOperator.reflectValue(t, -distances.backwardDistance, distances.forwardDistance);

        for (int i = 0; i < values.length; i++) {
            sample[i] = values[i] - sample[i] * t;
        }

        parameter.setAllParameterValuesQuietly(sample);
        parameter.fireParameterChangedEvent();

//        double tForward = t / (distances.forwardDistance * window);
//        double forwardLogDensity = uniformProductLogPdf(tForward);
//
//        double backWardDistance = distances.backwardDistance + t;
//        double tBackward = t / (backWardDistance * window);
//        double backwardLogDensity = uniformProductLogPdf(tBackward);

        return 0;
    }

//    private double uniformProductLogPdf(double t) {
//        double density = -Math.log(t);
//        return Math.log(density);
//    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (value > 0) value = 0;
        window = Math.exp(value);
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(window);
    }

    @Override
    public double getRawParameter() {
        return window;
    }

    @Override
    public String getAdaptableParameterName() {
        return WINDOW_SIZE;
    }

    @Override
    public String getOperatorName() {
        return CONVEX_RW;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }


}
