package dr.inference.operators;

import dr.inference.distribution.*;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.GammaDistribution;
import dr.util.Attribute;
import dr.xml.*;

/**
 * Created by mkarcher on 4/12/17.
 */
public class RandomWalkGammaPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    public static final String OPERATOR_NAME = "randomWalkGammaPrecisionGibbsOperator";
    public static final String DATA = "data";
    public static final String PRECISION = "precision";
    public static final String PRIOR = "prior";

    public RandomWalkGammaPrecisionGibbsOperator(Parameter data, Parameter precision, Distribution prior, double weight) {
        if (!(prior instanceof GammaDistribution || prior instanceof GammaDistributionModel))
            throw new RuntimeException("Precision prior must be Gamma");

        this.precision = precision;
        this.data = data;
        this.prior = prior;

        setWeight(weight);
    }

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     * @throws OperatorFailedException if operator fails and should be rejected
     */
    public double doOperation() throws OperatorFailedException {

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

        // Calculate weighted sum-of-squares
//        final double mu = meanParameter.getParameterValue(0);
        double SSE = 0;
        double x;
        int n = data.getDimension();
        for (int i = 1; i < n; i++) {
            x = data.getParameterValue(i) - data.getParameterValue(i-1);
            SSE += x * x;
        }

        final double shape = priorShape + n / 2.0;
        final double rate = priorRate + 0.5 * SSE;

        final double draw = MathUtils.nextGamma(shape, rate); // Gamma( \alpha + n/2 , \beta + (1/2)*SSE )
        precision.setParameterValue(0, draw);

        return 0;
    }

    /**
     * @return the number of steps the operator performs in one go.
     */
    public int getStepCount() {
        return 1;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            Parameter data = (Parameter) ((XMLObject) xo.getChild(DATA)).getChild(Parameter.class);
            Parameter precision = (Parameter) ((XMLObject) xo.getChild(PRECISION)).getChild(Parameter.class);
            DistributionLikelihood prior = (DistributionLikelihood) ((XMLObject) xo.getChild(PRIOR)).getChild(DistributionLikelihood.class);

//            System.err.println("class: " + prior.getDistribution().getClass());

            if (!((prior.getDistribution() instanceof GammaDistribution) ||
                    (prior.getDistribution() instanceof GammaDistributionModel)) )
                throw new XMLParseException("Gibbs operator assumes normal-gamma model");

            return new RandomWalkGammaPrecisionGibbsOperator(data, precision, prior.getDistribution(), weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an operator on the precision parameter of a normal random walk model with gamma prior.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
                new ElementRule(PRECISION,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
                new ElementRule(PRIOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
        };

    };

    private final Distribution prior;
    private final Parameter data;
    private final Parameter precision;
}
