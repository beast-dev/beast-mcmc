package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.distributions.ConvexSpaceRandomGenerator;
import jebl.math.Random;

public class ConvexSpaceRandomWalkOperator extends AbstractAdaptableOperator {
    private double window;
    private final ConvexSpaceRandomGenerator generator;
    private final Parameter parameter;

    public static final String CONVEX_RW = "convexSpaceRandomWalkOperator";
    public static final String WINDOW_SIZE = "relativeWindowSize";

    public ConvexSpaceRandomWalkOperator(Parameter parameter, ConvexSpaceRandomGenerator generator,
                                         double window, double weight) {
        setWeight(weight);

        this.parameter = parameter;
        this.generator = generator;
        this.window = window;
    }


    @Override
    public double doOperation() {
        double[] sample = (double[]) generator.nextRandom();
        double[] values = parameter.getParameterValues();

        ConvexSpaceRandomGenerator.LineThroughPoints distances = generator.distanceToEdge(values, sample);

        double t = window * Random.nextDouble();
        double oneMinus = 1.0 - t;

        for (int i = 0; i < values.length; i++) {
            sample[i] = values[i] * oneMinus + sample[i] * t;
        }

        parameter.setAllParameterValuesQuietly(sample);
        parameter.fireParameterChangedEvent();

        double tForward = t / (distances.forwardDistance * window);
        double forwardLogDensity = uniformProductLogPdf(tForward);

        double backWardDistance = distances.backwardDistance + t;
        double tBackward = t / (backWardDistance * window);
        double backwardLogDensity = uniformProductLogPdf(tBackward);

        return backwardLogDensity - forwardLogDensity;
    }

    private double uniformProductLogPdf(double t) {
        double density = -Math.log(t);
        return Math.log(density);
    }


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
