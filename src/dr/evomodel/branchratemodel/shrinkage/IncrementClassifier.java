package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.inference.distribution.ParametricMultivariateDistributionModel;

import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgePriorSampler;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider, Loggable {
    //todo: add logger
    public static final String INCREMENT_CLASSIFIER = "incrementClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";
//    public static final String SAMPLE_PRIOR = "samplePrior";

    private AutoCorrelatedBranchRatesDistribution acbr; //autocorrelated branch rates
    private DifferentiableBranchRates branchRateModel;
    private double epsilon;
    //    private double targetProb;
    private int dim;
    private Helper helper;

    private double sd;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution acbr, double epsilon) {
        this.acbr = acbr;
        this.branchRateModel = acbr.getBranchRateModel();
        this.epsilon = epsilon;
//        this.sd = sd; //only gets used if we are using targetProb

//        if (epsilon != 0.0) {
//            this.epsilon = epsilon;
//        } else if (targetProb != 0.0) {
//            this.targetProb = targetProb;
//            this.epsilon = findEpsilonFromTargetProb(targetProb);
//        }

        this.dim = acbr.getDimension();
        this.classified = new double[dim];
        classify();

        helper = new Helper();
        setupTraits();
    }

//    private double findEpsilonFromTargetProb(double targetProb) {
//        // finds eps such that Prob{abs(increment) < eps} = targetProb
//        double probToFindEpsilon = (targetProb / 2) + 0.5;
////        NormalDistribution normal = new NormalDistribution(0, sd);
////        return normal.quantile(probToFindEpsilon);
//        return()
//    }

    private void classify() {
        double increment;
        for (int i = 0; i < dim; i++) {
            increment = acbr.getIncrement(i);
            classified[i] = 1; // is an increment
            if (Math.abs(increment) < epsilon) {
                classified[i] = 0; // is not an increment
            }
        }
    }

    private void setupTraits() {
        TreeTrait classifierTrait = new TreeTrait.D() {

            public String getTraitName() {
                return "incrementClassifier";
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public Double getTrait(Tree tree, NodeRef node) {
                classify();
                int index = branchRateModel.getParameterIndexFromNode(node);
                return classified[index];
            }

            public boolean getLoggable() {
                return true;
            }
        };
        helper.addTrait(classifierTrait);
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INCREMENT_CLASSIFIER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            AutoCorrelatedBranchRatesDistribution acbr = (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);

            if (xo.hasAttribute(EPSILON) && xo.hasAttribute(TARGET_PROBABILITY)) {
                throw new XMLParseException("Cannot set both epsilon and target probability.");
            }

//            boolean samplePrior = xo.getAttribute(SAMPLE_PRIOR, false);

//            double sd = 1.0; //todo: put magic number here.
//            if (xo.hasAttribute(TARGET_PROBABILITY) && samplePrior){
//                System.out.println("Using default Bayesian bridge prior with standard deviation " + sd + " to set epsilon.");
//            }

            double epsilon = xo.getAttribute(EPSILON, 0.0);
            double targetProbability = xo.getAttribute(TARGET_PROBABILITY, 0.0);


            BayesianBridgePriorSampler bayesianBridgePriorSampler;
            if (targetProbability != 0.0) {
                bayesianBridgePriorSampler = (BayesianBridgePriorSampler) xo.getChild(BayesianBridgePriorSampler.class);
//                sd = bayesianBridgePriorSampler.getStandardDeviation();
                epsilon = bayesianBridgePriorSampler.getEpsilon(targetProbability);
            }

            if (epsilon < 0.0) {
                throw new XMLParseException("epsilon must be positive.");
            } else if (targetProbability < 0.0) {
                throw new XMLParseException("target probability must be positive.");
            }

            IncrementClassifier incrementClassifier = new IncrementClassifier(acbr, epsilon);
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
            return "Classifies increment as 0 or 1 based on arbitrary density cutoff epsilon.";
        }

        @Override
        public Class getReturnType() {
            return IncrementClassifier.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(AutoCorrelatedBranchRatesDistribution.class),
//                new ElementRule(BayesianBridgePriorSampler.class, true),
//                AttributeRule.newDoubleRule(EPSILON, true),
//                AttributeRule.newDoubleRule(TARGET_PROBABILITY, true),
//                AttributeRule.newBooleanRule(SAMPLE_PRIOR, true),
                new XORRule(
                        AttributeRule.newDoubleRule(EPSILON, false),
                        new AndRule(
                                AttributeRule.newDoubleRule(TARGET_PROBABILITY, false),
                                new ElementRule(BayesianBridgePriorSampler.class, false)
                        )
                )
        };
    };

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }
}
