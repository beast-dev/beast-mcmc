package dr.evomodel.coalescent;

import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author joseph
 *         Date: 25/04/2008
 */

// It should be a model since it may be the only user of parameter sigma

public class BMPriorLikelihood extends AbstractModelLikelihood {
    // private final Parameter mean;
    private final Parameter sigma;
    //private final Parameter lambda;
    private final boolean logSpace;
    //private final boolean normalize;
    // private Parameter data;
    // private Parameter times;
    private ParametricDistributionModel popMeanPrior = null;
    private final VariableDemographicModel m;

    private BMPriorLikelihood(Parameter sigma, boolean logSpace, VariableDemographicModel m) {
        super(BM);

        this.logSpace = logSpace;
        //this.normalize = normalize;

        //this.mean = mean;
//        if( mean != null ) {
//            addParameter(mean);
//        }

        this.m = m;
        this.sigma = sigma;
        addParameter(sigma);

        // this.lambda = lambda;
//        if (lambda != null) {
//            addParameter( lambda );
//        }
    }

    public BMPriorLikelihood(Parameter sigma, VariableDemographicModel demographicModel, boolean logSpace,
                             ParametricDistributionModel popMeanPrior) {
        this(sigma, logSpace, demographicModel);
        //   this.m = demographicModel;
        addModel(demographicModel);

        //   this.data = this.times = null;
        this.popMeanPrior = popMeanPrior;
    }

//    public BMPriorLikelihood(Parameter sigma,
//                                            Parameter dataParameter, Parameter timesParameter, boolean logSpace, boolean normalize) {
//        this(sigma, logSpace);
//        dataParameter.addParameterListener(this);
//        timesParameter.addParameterListener(this);
//        this.data = dataParameter;
//        this.times = timesParameter;
//    }

    // log of normal distribution coeffcient.
    final private double logNormalCoef = -0.5 * Math.log(2 * Math.PI);

    // A specialized version where everything is normalized. Time is normalized to 1. Data moved to mean zero and rescaled
    // according to time.

    private double reNormalize(VariableDemographicModel m) {

        final double[] tps = m.getDemographicFunction().allTimePoints();
        // get a copy since we re-scale data
        final double[] vals = m.getPopulationValues().getParameterValues();

        //assert ! logSpace : "not implemented yet";

        final double tMax = tps[tps.length - 1];


        if (false) {
            return -Math.log(tMax / m.getDemographicFunction().getIntegral(0, tMax));
        } else {
            // compute mean
            double popMean = tMax / m.getDemographicFunction().getIntegral(0, tMax);
            // todo not correct when using midpoints
//        if( m.getType() == VariableDemographicModel.Type.LINEAR ) {
//            for(int k = 0; k < tps.length; ++k) {
//                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
//                popMean += dt * (vals[k+1] + vals[k]);
//            }
//            popMean /= (2* tMax);
//        } else {
//            for(int k = 0; k < tps.length; ++k) {
//                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
//                popMean += dt * vals[k];
//            }
//            popMean /= tMax;
//        }

            // Normalize to time interval = 1 and mean = 0
            final double sigma = this.sigma.getStatisticValue(0);
//        final double lam = 0.5 * tMax;
//        for(int k = 0; k < vals.length; ++k) {
//            vals[k] = vals[k]/tMax;
//        }

            // optimized version of the code in getLogLikelihood.
            // get factors out when possible. logpdf of a normal is -x^2/2, when mean is 0
            double ll = 0.0;

            final double s2 = 2 * sigma * sigma;
            for (int k = 0; k < tps.length; ++k) {
                final double dt = ((tps[k] - (k > 0 ? tps[k - 1] : 0)) / tMax);

                //final double d = (vals[k+1] - vals[k]);
                final double r = logSpace ? (vals[k + 1] - vals[k]) : Math.log(vals[k + 1] / vals[k]);
                final double d = r;
                ll -= (d * d / (s2 * dt));
                ll -= 0.5 * Math.log(dt);
            }
            ll += tps.length * (logNormalCoef - Math.log(sigma));
            //ll /= tps.length;
            if (popMeanPrior != null) {
                ll += popMeanPrior.logPdf(popMean);
            } else {
                // default Jeffrey's
                ll -= logSpace ? popMean : Math.log(popMean);
            }

            return ll;
        }
    }

    public double getLogLikelihood() {
//        if( lastValue != Double.NEGATIVE_INFINITY ) {
//            return lastValue;
//        }

        assert m != null;
        lastValue = reNormalize(m);
        return lastValue;
    }

    public void makeDirty() {

        lastValue = Double.NEGATIVE_INFINITY;
    }

    static final String BM = "BrownianMotion";
    static final String LOG_SPACE = "logUnits";
    //static final String NORMALIZE = "normalize";

    //static final String DATA = "data";
    // static final String TIMES = "times";

    // public static final String MEAN = "mean";
    public static final String SIGMA = "sigma";
    //  public static final String LAMBDA = "lambda";

    // (todo) Parser is still in a bad state
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return BMPriorLikelihood.class;
        }

        public String getParserName() {
            return BM;
        }

        private Parameter getParam(XMLObject xo, String name) throws XMLParseException {
            final XMLObject object = xo.getChild(name);
            // optional
            if (object == null) {
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

            // Parameter mean = getParam(xo, MEAN);
            Parameter sigma = getParam(xo, SIGMA);
            //  Parameter lambda = getParam(xo, LAMBDA);

            final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
            //   final boolean normalize = xo.getAttribute(NORMALIZE, false);

            VariableDemographicModel m = (VariableDemographicModel) xo.getChild(VariableDemographicModel.class);

            ParametricDistributionModel popMeanPrior = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);
            return new BMPriorLikelihood(sigma, m, logSpace, popMeanPrior);
        }

        public XMLSyntaxRule[] getSyntaxRules() {


            return new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(LOG_SPACE, true),

                    new XORRule(
                            new ElementRule(SIGMA, Double.class),
                            new ElementRule(SIGMA, Parameter.class)
                    ),

                    // you can't have a XOR (b AND c), yikes
                    // make all optional and check in parser
                    new ElementRule(VariableDemographicModel.class),

//                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true),
//                    new ElementRule(TIMES, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true)

            };
        }
    };

    // simply saves last value
    double lastValue = Double.NEGATIVE_INFINITY;

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
        makeDirty();
    }

    protected void acceptState() {
    }
}