package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.distributions.RandomGenerator;
import jebl.math.Random;

public class ConvexSpaceRandomWalkOperator extends AbstractAdaptableOperator {
    private double window;
    private final RandomGenerator generator;
    private final Parameter parameter;

    public static final String CONVEX_RW = "convexSpaceRandomWalkOperator";
    public static final String WINDOW_SIZE = "relativeWindowSize";

    public ConvexSpaceRandomWalkOperator(Parameter parameter, RandomGenerator generator, double window, double weight) {
        setWeight(weight);

        this.parameter = parameter;
        this.generator = generator;
        this.window = window;
    }


    @Override
    public double doOperation() {
        double[] sample = (double[]) generator.nextRandom();
        double[] values = parameter.getParameterValues();
        double t = window * Random.nextDouble();
        double oneMinus = 1.0 - t;

        for (int i = 0; i < values.length; i++) {
            sample[i] = values[i] * t + sample[i] * oneMinus;
        }

        parameter.setAllParameterValuesQuietly(sample);
        parameter.fireParameterChangedEvent();
        return 0.0; //TODO: need to check that the prior is uniform
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


}
