package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class NormalNormalMeanGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String OPERATOR_NAME = "normalNormalMeanGibbsOperator";
    public static final String LIKELIHOOD = "likelihood";
    public static final String PRIOR = "prior";

    public NormalNormalMeanGibbsOperator(DistributionLikelihood inLikelihood, Distribution prior,
                                         double weight) {

        if (!(prior instanceof NormalDistribution || prior instanceof NormalDistributionModel))
            throw new RuntimeException("Mean prior must be Normal");

        this.likelihood = inLikelihood.getDistribution();
        this.dataList = inLikelihood.getDataList();
        if (likelihood instanceof NormalDistributionModel)
            this.meanParameter = ((NormalDistributionModel) likelihood).getMeanParameter();
        else if (likelihood instanceof LogNormalDistributionModel) {
            this.meanParameter = ((LogNormalDistributionModel) likelihood).getMeanParameter();
            isLog = true;
        } else
            throw new RuntimeException("Likelihood must be Normal or log Normal");

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

        double priorPrecision = 1.0 / prior.variance();
        double priorMean = prior.mean();
        double likelihoodPrecision = 1.0 / likelihood.variance();

        double total = 0;
        int n = 0;
        for (Statistic statistic : dataList) {
            for (double x : statistic.getAttributeValue()) {
                if (isLog)
                    total += Math.log(x);
                else
                    total += x;
                n++;
            }
        }

        double precision = priorPrecision + likelihoodPrecision * n;
        double mu = (priorPrecision * priorMean + likelihoodPrecision * total) / precision;
        meanParameter.setParameterValue(0,
                MathUtils.nextGaussian() / Math.sqrt(precision) + mu);  // N(\mu, \precision)
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

            DistributionLikelihood likelihood = (DistributionLikelihood) ((XMLObject) xo.getChild(LIKELIHOOD)).getChild(DistributionLikelihood.class);
            DistributionLikelihood prior = (DistributionLikelihood) ((XMLObject) xo.getChild(PRIOR)).getChild(DistributionLikelihood.class);

            System.err.println("class: " + prior.getDistribution().getClass());

            if (!((prior.getDistribution() instanceof NormalDistribution) ||
                    (prior.getDistribution() instanceof NormalDistributionModel)
            ) ||
                    !((likelihood.getDistribution() instanceof NormalDistributionModel) ||
                            (likelihood.getDistribution() instanceof LogNormalDistributionModel)
                    ))
                throw new XMLParseException("Gibbs operator assumes normal-normal model");

            return new NormalNormalMeanGibbsOperator(likelihood, prior.getDistribution(), weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a operator on the mean parameter of a normal model with normal prior.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(LIKELIHOOD,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
                new ElementRule(PRIOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
        };

    };

    private Distribution likelihood;
    private Distribution prior;
    private boolean isLog = false;

    private List<Statistic> dataList;
    private Parameter meanParameter;
}