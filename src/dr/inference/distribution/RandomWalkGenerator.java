package dr.inference.distribution;

import dr.inference.model.*;
import dr.inferencexml.distribution.RandomWalkGeneratorParser;
import dr.math.MathUtils;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.NormalDistribution;

/**
 * Created by mkarcher on 4/3/17.
 */
public class RandomWalkGenerator extends AbstractModelLikelihood implements GaussianProcessRandomGenerator {

    public RandomWalkGenerator(Parameter data, Parameter firstElementPrecision, Parameter precision) {
        super(RandomWalkGeneratorParser.RANDOM_WALK_GENERATOR);
        this.data = data;
        this.dimension = data.getDimension();
        this.firstElementPrecision = firstElementPrecision;
        this.precision = precision;
//        this.logScale = logScale;
    }

    @Override
    public double[] nextRandom() {
        double[] result = new double[dimension];
        result[0] = MathUtils.nextGaussian() / Math.sqrt(firstElementPrecision.getParameterValue(0));
        for (int i = 1; i < dimension; i++) {
            double eps = MathUtils.nextGaussian();
            result[i] = result[i-1] + eps / Math.sqrt(precision.getParameterValue(0));
        }
        return result;
    }

    @Override
    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }

    public double logPdf(double[] x) {
        double result = 0;

        result += NormalDistribution.logPdf(x[0], 0, 1/Math.sqrt(firstElementPrecision.getParameterValue(0)));
        for (int i = 1; i < dimension; i++) {
            result += NormalDistribution.logPdf(x[i], x[i-1], 1/Math.sqrt(precision.getParameterValue(0)));
        }

        return result;
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        return new double[0][];
    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        return logPdf(data.getParameterValues());
    }

    @Override
    public void makeDirty() {

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    private final Parameter data;
    private int dimension;
    private Parameter firstElementPrecision;
    private Parameter precision;
//    private final boolean logScale;
}
