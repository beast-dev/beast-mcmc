package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.ParametricMultivariateDistributionModel;

import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider {
// todo: find magic default s.d. from prior simulation (then use targetProb to find eps)
    public static final String PARAMETER_CLASSIFIER = "parameterClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";

    private AutoCorrelatedBranchRatesDistribution brm;
    private double epsilon;
    private boolean epsilonKnown = false; // todo: cleanup -- probably don't need to dynamically update epsilon
    private double targetProb;
    private int dim;

    private double sd = 1.0; //todo: put magic number here.
//    private ParametricMultivariateDistributionModel distribution;
    private JointBayesianBridgeDistributionModel distribution;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution brm, double epsilon, double targetProb) {
        this.brm = brm;
        this.distribution = (distribution instanceof JointBayesianBridgeDistributionModel) ?
                (JointBayesianBridgeDistributionModel) distribution : null;

        if (epsilon != 0.0) {
            this.epsilon = epsilon;
            epsilonKnown = true;
        } else if (targetProb != 0.0) {
            this.targetProb = targetProb;
            this.epsilon = findEpsilonFromTargetProb(targetProb);
        }

        this.dim = brm.getDimension();
        this.classified = new double[dim];
        classify();
    }

    private double findEpsilonFromTargetProb(double targetProb) {
        // finds eps such that Prob{abs(increment) < eps} = targetProb
        double probToFindEpsilon = (targetProb / 2) + 0.5;
        NormalDistribution norm = new NormalDistribution(0, sd);
        epsilonKnown = true;
        return norm.quantile(probToFindEpsilon);
    }

    private void classify() {
        double increment;
        if(!epsilonKnown){
            findEpsilonFromTargetProb(targetProb);
        }
        for (int i = 0; i < dim; i++) {
            increment = brm.getIncrement(i);
            classified[i] = 1; // is an increment
            if (Math.abs(increment) < epsilon) {
                classified[i] = 0; // is not an increment
            }
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[0];
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return null;
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARAMETER_CLASSIFIER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            AutoCorrelatedBranchRatesDistribution brm = (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);

            if (xo.hasAttribute(EPSILON) && xo.hasAttribute(TARGET_PROBABILITY)) {
                throw new XMLParseException("Cannot set both epsilon and target probability.");
            }

            double epsilon = xo.getAttribute(EPSILON, 0);
            double targetProbability = xo.getAttribute(TARGET_PROBABILITY, 0);

            if (epsilon < 0) {
                throw new XMLParseException("epsilon must be positive.");
            }
            else if (targetProbability < 0) {
                throw new XMLParseException("target probability must be positive.");
            }

            IncrementClassifier incrementClassifier = new IncrementClassifier(brm, epsilon, targetProbability);
            return incrementClassifier;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Reports 0 or 1 depending on whether a model parameter falls within a specified density range.";
        }

        @Override
        public Class getReturnType() {
            return IncrementClassifier.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                AttributeRule.newDoubleRule(EPSILON, true),
                AttributeRule.newDoubleRule(TARGET_PROBABILITY, true),
        };
    };
}
