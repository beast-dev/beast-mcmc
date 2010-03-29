package dr.inference.operators;

import dr.inference.distribution.*;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.Distribution;
import dr.math.MathUtils;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class RegressionGibbsPrecisionOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "regressionGibbsPrecisionOperator";

    private LinearRegression linearModel;
    private Parameter precision;
    private int dim;
    private int N;
    private int[] scaleDesign;
    private Distribution prior;

    public RegressionGibbsPrecisionOperator(LinearRegression linearModel, Parameter precision, Distribution prior) {

        super();
        if (!(prior instanceof GammaDistribution || prior instanceof GammaDistributionModel))
              throw new RuntimeException("Precision prior must be Gamma");
        this.prior = prior;
        this.linearModel = linearModel;
        this.precision = precision;
        this.dim = precision.getDimension();
        scaleDesign = linearModel.getScaleDesign();
        N = linearModel.getDependentVariable().getDimension();
    }

    public int getStepCount() {
        return 1;
    }

    public double doOperation() throws OperatorFailedException {

        double[] Y = linearModel.getTransformedDependentParameter();
        double[] xBeta = linearModel.getXBeta();

        final double priorMean = prior.mean();
        final double priorVariance = prior.variance();

        double priorRate;
        double priorShape;

        if (priorMean == 0) {
            priorRate = 0;
            priorShape = -0.5; // Uninformative prior
        } else {
            priorRate = priorMean / priorVariance;
            priorShape = priorMean * priorRate;
        }

        for (int k = 0; k < dim; k++) { // Do draw for precision[k]

            // Calculate weighted sum-of-squares
            double SSE = 0;
            int n = 0;

            for(int i=0; i<N; i++) {
                if(scaleDesign[i] == k) {
                    SSE += (Y[i] - xBeta[i])*(Y[i] - xBeta[i]);
                    n++;
                }
            }

            final double shape = priorShape + n / 2.0;
            final double rate = priorRate + 0.5 * SSE;

            final double draw = MathUtils.nextGamma(shape, rate); // Gamma( \alpha + n/2 , \beta + (1/2)*SSE )
            precision.setParameterValue(k, draw);
        }

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            LinearRegression linearModel = (LinearRegression) xo.getChild(LinearRegression.class);
            Parameter precision = (Parameter) xo.getChild(Parameter.class);
            DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);

            if (!((prior.getDistribution() instanceof GammaDistribution) ||
                        (prior.getDistribution() instanceof GammaDistributionModel)
                ))
                    throw new XMLParseException("Gibbs operator assumes normal-gamma model");

            RegressionGibbsPrecisionOperator operator = new RegressionGibbsPrecisionOperator(linearModel,precision,prior.getDistribution());
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on an internal node trait.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(Parameter.class),
                new ElementRule(DistributionLikelihood.class),
                new ElementRule(LinearRegression.class),
        };

    };

}