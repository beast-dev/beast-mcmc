package dr.math.distributions;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.*;

public class MultivariateGammaLikelihood extends AbstractModelLikelihood implements ParametricMultivariateDistributionModel, GradientWrtParameterProvider {

    private final Parameter shape;
    private final Parameter scale;
    private final Parameter data;
    private final int dim;
    public static final String TYPE = "multivariateGamma";


    public MultivariateGammaLikelihood(String name, Parameter shape, Parameter scale, Parameter data) {
        super(name);
        this.data = data;
        this.shape = shape;
        this.scale = scale;
        this.dim = data.getDimension();
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("not implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // do nothing
    }

    @Override
    protected void storeState() {
        // do nothing
    }

    @Override
    protected void restoreState() {
        // do nothing
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    @Override
    public double logPdf(double[] x) {
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            sum += GammaDistribution.logPdf(x[i], shape.getParameterValue(i), scale.getParameterValue(i));
        }
        return sum;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("not yet implemented");

    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return data;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] grad = new double[dim];
        for (int i = 0; i < dim; i++) {
            grad[i] = GammaDistribution.gradLogPdf(data.getParameterValue(i), shape.getParameterValue(i),
                    scale.getParameterValue(i));
        }
        return grad;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double[] values = data.getParameterValues();
        return logPdf(values);
    }

    @Override
    public void makeDirty() {
        // do nothing
    }

    private static final String SCALE = "scale";
    private static final String SHAPE = "shape";
    private static final String DATA = "data";
    private static final String MULTIVARIATE_GAMMA_LIKELIHOOD = "multivariateGammaLikelihood";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);
            Parameter shape = (Parameter) xo.getChild(SHAPE).getChild(Parameter.class);
            Parameter scale = (Parameter) xo.getChild(SCALE).getChild(Parameter.class);

            if (data.getDimension() != shape.getDimension() || data.getDimension() != scale.getDimension()) {
                throw new XMLParseException("All parameters must have the same dimension. " +
                        "The components currently have the following dimensions:\n" +
                        "\t" + DATA + ": " + data.getDimension() + "\n" +
                        "\t" + SHAPE + ": " + shape.getDimension() + "\n" +
                        "\t" + SCALE + ": " + scale.getDimension() + "\n");
            }

            return new MultivariateGammaLikelihood(MULTIVARIATE_GAMMA_LIKELIHOOD, shape, scale, data);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(SCALE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(SHAPE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            };
        }

        @Override
        public String getParserDescription() {
            return "independent gamma distributions with common parameter";
        }

        @Override
        public Class getReturnType() {
            return MultivariateGammaLikelihood.class;
        }

        @Override
        public String getParserName() {
            return MULTIVARIATE_GAMMA_LIKELIHOOD;
        }
    };
}
