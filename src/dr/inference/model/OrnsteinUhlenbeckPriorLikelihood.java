package dr.inference.model;

import dr.evomodel.coalescent.VariableDemographicModel;
import dr.math.NormalDistribution;
import dr.xml.*;

/**
 * Ornstein-Uhlenbeck prior.
 * <p/>
 * A diffusion process - Basically a correlated sequence of normally distributed values, where correlation is time
 * dependent.
 * <p/>
 * <p/>
 * Very experimental and only slightly tested at this time
 *
 * @author joseph
 *         Date: 25/04/2008
 */
public class OrnsteinUhlenbeckPriorLikelihood extends Likelihood.Abstract implements ParameterListener {
    private Parameter mean;
    private Parameter sigma;
    private Parameter lambda;
    private boolean logSpace;
    private boolean normalize;
    private Parameter data;
    private Parameter times;

    private OrnsteinUhlenbeckPriorLikelihood(Model m, Parameter mean, Parameter sigma, Parameter lambda,
                                             boolean logSpace, boolean normalize) {
        super(m);
        mean.addParameterListener(this);
        sigma.addParameterListener(this);

        this.logSpace = logSpace;
        this.normalize = normalize;

        this.mean = mean;
        this.sigma = sigma;
        if (lambda != null) {
            this.lambda = lambda;
            lambda.addParameterListener(this);
        } else {
            this.lambda = new Parameter.Default(0.5);
        }
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            VariableDemographicModel demographicModel, boolean logSpace, boolean normalize) {
        this(demographicModel, mean, sigma, lambda, logSpace, normalize);
        this.data = this.times = null;
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            Parameter dataParameter, Parameter timesParameter, boolean logSpace, boolean normalize) {
        this(null, mean, sigma, lambda, logSpace, normalize);
        dataParameter.addParameterListener(this);
        timesParameter.addParameterListener(this);
        this.data = dataParameter;
        this.times = timesParameter;
    }

    protected double calculateLogLikelihood() {

        VariableDemographicModel m = (VariableDemographicModel) getModel();
        final double[] tps = m != null ? m.getDemographicFunction().allTimePoints() : times.getParameterValues();

        final double[] vals = m != null ? m.getPopulationValues().getParameterValues() : data.getParameterValues();

        if (logSpace) {
            for (int k = 0; k < vals.length; ++k) {
                vals[k] = Math.log(vals[k]);
            }
        }

        final double lambda = this.lambda.getStatisticValue(0);
        final double mean = this.mean.getStatisticValue(0);

        double sigma = this.sigma.getStatisticValue(0);
        if (normalize) {
            // make the process have a SD of sigma
            sigma *= Math.sqrt(2 * lambda);
        }


        double logL = NormalDistribution.logPdf(vals[0], mean, sigma);

        final double xScale = -lambda * (normalize ? 1.0 / tps[tps.length - 1] : 1.0);

        for (int k = 0; k < tps.length; ++k) {
            double dt = tps[k] - (k > 0 ? tps[k - 1] : 0);
            double a = Math.exp(xScale * dt);
            double den = sigma * Math.sqrt((1 - a * a) / (2 * lambda));
            double z = (vals[k + 1] - (vals[k] * a + mean * (1 - a))) / den;
            logL += NormalDistribution.logPdf(z, 0, 1);
        }
        return logL;
    }

    public void parameterChangedEvent(Parameter parameter, int index) {
        makeDirty();
    }

    static final String OU = "Ornstein-Uhlenbeck";
    static final String LOG_SPACE = "logUnits";
    static final String NORMALIZE = "normalize";

    static final String DATA = "data";
    static final String TIMES = "times";

    public static final String MEAN = "mean";
    public static final String SIGMA = "sigma";
    public static final String LAMBDA = "lambda";


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return OrnsteinUhlenbeckPriorLikelihood.class;
        }

        public String getParserName() {
            return OU;
        }

        private Parameter getParam(XMLObject xo, String name) throws XMLParseException {
            final XMLObject object = (XMLObject) xo.getChild(name);
            final Object child = object.getChild(0);
            if (child instanceof Parameter) {
                return (Parameter) child;
            }

            double x = object.getDoubleChild(0);
            return new Parameter.Default(x);
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter mean = getParam(xo, MEAN);
            Parameter sigma = getParam(xo, SIGMA);
            Parameter lambda = getParam(xo, LAMBDA);

            final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
            final boolean normalize = xo.getAttribute(NORMALIZE, false);

            VariableDemographicModel m = (VariableDemographicModel) xo.getChild(VariableDemographicModel.class);
            if (m != null) {
                return new OrnsteinUhlenbeckPriorLikelihood(mean, sigma, lambda, m, logSpace, normalize);
            }

            XMLObject cxo1 = (XMLObject) xo.getChild(DATA);
            Parameter dataParameter = (Parameter) cxo1.getChild(Parameter.class);
            XMLObject cxo2 = (XMLObject) xo.getChild(TIMES);
            Parameter timesParameter = (Parameter) cxo2.getChild(Parameter.class);
            return new OrnsteinUhlenbeckPriorLikelihood(mean, sigma, lambda, dataParameter, timesParameter, logSpace, normalize);
        }

        public XMLSyntaxRule[] getSyntaxRules() {


            return new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(LOG_SPACE, true),
                    AttributeRule.newBooleanRule(NORMALIZE, true),
                    //new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
                    //new ElementRule(TIMES, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
                    new XORRule(
                            new ElementRule(MEAN, Double.class),
                            new ElementRule(MEAN, Parameter.class)
                    ),
                    new XORRule(
                            new ElementRule(SIGMA, Double.class),
                            new ElementRule(SIGMA, Parameter.class)
                    ),

                    new XORRule(
                            new ElementRule(LAMBDA, Double.class),
                            new ElementRule(LAMBDA, Parameter.class)
                    ),

                    // you can't have a XOR (b AND c), yikes
                    // make all optional and check in parser
                    new ElementRule(VariableDemographicModel.class, true),

                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true),
                    new ElementRule(TIMES, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true)

            };
        }

    };
}