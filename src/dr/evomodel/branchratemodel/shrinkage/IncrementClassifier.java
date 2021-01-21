package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider {
// todo: find magic default epsilon from prior simulation
    public static final String PARAMETER_CLASSIFIER = "parameterClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";

    private AutoCorrelatedBranchRatesDistribution brm;
    private double epsilon;
    private boolean epsilonKnown = false;
    private double targetProb;
    private int dim;
    private ParametricMultivariateDistributionModel distribution;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution brm, double epsilon, double targetProb) {
        this.brm = brm;
        this.distribution = brm.getPrior();

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
        // finds eps such that p(abs(increment)) < epsilon = targetProb
        throw new RuntimeException("Target probabiltiy not yet implemented.");
    }

    private void classify() {
        double increment;
        if(!epsilonKnown){
            findEpsilonFromTargetProb(targetProb);
        }
        for (int i = 0; i < dim; i++) {
            increment = brm.getIncrement(i);
            classified[i] = 1; // is an increment
            if (increment < epsilon) {
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
