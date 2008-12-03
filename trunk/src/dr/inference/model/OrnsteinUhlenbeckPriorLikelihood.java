package dr.inference.model;

import dr.evomodel.coalescent.VariableDemographicModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.math.distributions.NormalDistribution;
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

// It should be a model since it may be the only user of parameter sigma

public class OrnsteinUhlenbeckPriorLikelihood extends AbstractModel implements Likelihood {
    private final Parameter mean;
    private final Parameter sigma;
    private final Parameter lambda;
    private final boolean logSpace;
    private final boolean normalize;
    private Parameter data;
    private Parameter times;
    private ParametricDistributionModel popMeanPrior = null;
    private VariableDemographicModel m = null;

    private OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                             boolean logSpace, boolean normalize) {
        super(OU);

        this.logSpace = logSpace;
        this.normalize = normalize;

        this.mean = mean;
        if( mean != null ) {
            addParameter(mean);
        }

        this.sigma = sigma;
        addParameter(sigma);

        this.lambda = lambda;
        if (lambda != null) {
            addParameter( lambda );
        }
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            VariableDemographicModel demographicModel, boolean logSpace, boolean normalize,
                                            ParametricDistributionModel popMeanPrior) {
        this(mean, sigma, lambda, logSpace, normalize);
        this.m = demographicModel;

        this.data = this.times = null;
        this.popMeanPrior = popMeanPrior;
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            Parameter dataParameter, Parameter timesParameter, boolean logSpace, boolean normalize) {
        this(mean, sigma, lambda, logSpace, normalize);
        dataParameter.addParameterListener(this);
        timesParameter.addParameterListener(this);
        this.data = dataParameter;
        this.times = timesParameter;
    }

    // log of normal distribution coeffcient.
    final private double logNormalCoef =  -0.5 * Math.log(2*Math.PI);

    // A specialized version where everthing is normalized. Time is normalized to 1. Data moved to mean zero and rescaled
    // according to time. Lambda is 0.5. The prior on mean is added.

    private double reNormalize(VariableDemographicModel m) {
        final double[] tps = m.getDemographicFunction().allTimePoints();
        // get a copy since we re-scale data
        final double[] vals = m.getPopulationValues().getParameterValues();

        assert ! logSpace : "not implemented yet";

        final double len = tps[tps.length-1];

        // compute mean
        double popMean = 0;
        // todo not correct when using midpoints
        if( m.getType() == VariableDemographicModel.Type.LINEAR ) {
            for(int k = 0; k < tps.length; ++k) {
                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
                popMean += dt * (vals[k+1] + vals[k]);
            }
            popMean /= (2* len);
        } else {
            for(int k = 0; k < tps.length; ++k) {
                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
                popMean += dt * vals[k];
            }
            popMean /= len;
        }

        // Normalize to time interval = 1 and mean = 0
        final double sigma = this.sigma.getStatisticValue(0)/ Math.sqrt(len);
        final double lam = 0.5 * len;
        for(int k = 0; k < vals.length; ++k) {
            vals[k] = (vals[k] - popMean)/len;
        }

        // optimized version of the code in getLogLikelihood.
        // get factors out when possible. logpdf of a normal is -x^2/2, when mean is 0
        double ll = 0.0;

        for(int k = 0; k < tps.length; ++k) {
            final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0)) / len;

            final double a = Math.exp(-lam * dt);
            final double d = (vals[k+1] - vals[k] * a);
            ll +=  (d*d) / (1-a*a);
        }

        final double f2 =  (2*lam) / (sigma*sigma);
        ll = tps.length * logNormalCoef + ll*f2/-2;

        if( popMeanPrior != null ) {
            ll += popMeanPrior.logPdf(popMean);
        } else {
            // default Jeffreys
            ll -= Math.log(popMean);
        }
        return ll;
    }

    public double getLogLikelihood() {
        if( lastValue > 0 ) {
            return lastValue;
        }
         double logL;

        if( normalize ) {
            assert m != null;
            logL = reNormalize(m);
        }  else {

//        final double[] tps = m != null ? m.getDemographicFunction().allTimePoints() : times.getParameterValues();
//
//        final double[] vals = m != null ? m.getPopulationValues().getParameterValues() : data.getParameterValues();

            final double[] tps = times.getParameterValues();
            final double[] vals = data.getParameterValues();

            if (logSpace) {
                for (int k = 0; k < vals.length; ++k) {
                    vals[k] = Math.log(vals[k]);
                }
            }

            final double lambda = this.lambda.getStatisticValue(0);
            final double mean = this.mean.getStatisticValue(0);

            double sigma = this.sigma.getStatisticValue(0);
            if( normalize ) {
                // make the process have a SD of sigma
                sigma *= Math.sqrt(2 * lambda);
            }

            logL = NormalDistribution.logPdf(vals[0], mean, sigma);

            final double xScale = -lambda * (normalize ? 1.0 / tps[tps.length - 1] : 1.0);

            for(int k = 0; k < tps.length; ++k) {
                double dt = tps[k] - (k > 0 ? tps[k - 1] : 0);
                double a = Math.exp(xScale * dt);
                double den = sigma * Math.sqrt((1 - a * a) / (2 * lambda));
                double z = (vals[k + 1] - (vals[k] * a + mean * (1 - a))) / den;
                logL += NormalDistribution.logPdf(z, 0, 1);
            }
        }

        lastValue = logL;

        return logL;
    }

    public void makeDirty() {
        lastValue = -1;
    }

    static final String OU = "Ornstein-Uhlenbeck";
    static final String LOG_SPACE = "logUnits";
    static final String NORMALIZE = "normalize";

    static final String DATA = "data";
    static final String TIMES = "times";

    public static final String MEAN = "mean";
    public static final String SIGMA = "sigma";
    public static final String LAMBDA = "lambda";

    // (todo) Parser is still in a bad state
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
            // optional
            if( object == null ) {
                return null;
            }
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
            
            if ( m != null ) {
                ParametricDistributionModel popMeanPrior = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);
                return new OrnsteinUhlenbeckPriorLikelihood(mean, sigma, lambda, m, logSpace, normalize, popMeanPrior);
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
                            new ElementRule(MEAN, Parameter.class),
                            true
                    ),
                    new XORRule(
                            new ElementRule(SIGMA, Double.class),
                            new ElementRule(SIGMA, Parameter.class)
                    ),

                    new XORRule(
                            new ElementRule(LAMBDA, Double.class),
                            new ElementRule(LAMBDA, Parameter.class),
                            true
                    ),

                    // you can't have a XOR (b AND c), yikes
                    // make all optional and check in parser
                    new ElementRule(VariableDemographicModel.class, true),

                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true),
                    new ElementRule(TIMES, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true)

            };
        }

    };

    // simply saves last value
    double lastValue = -1;

    public Model getModel() {
        return this;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[] {
                new dr.inference.loggers.NumberColumn(getId()) {
                    public double getDoubleValue() {
                        return getLogLikelihood();
                    }
                }
        };
    }
}