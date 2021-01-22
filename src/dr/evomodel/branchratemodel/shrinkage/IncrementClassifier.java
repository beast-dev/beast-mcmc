package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.inference.distribution.ParametricMultivariateDistributionModel;

import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider {
    // todo: find magic number based on input exponent and slabwidth -- sample from global and local scale
// todo: find magic default s.d. from prior simulation (then use targetProb to find eps)
    public static final String INCREMENT_CLASSIFIER = "incrementClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";

    private AutoCorrelatedBranchRatesDistribution acbr; //autocorrelated branch rates
    private DifferentiableBranchRates branchRateModel;
    private double epsilon;
    private double targetProb;
    private int dim;
    private Helper helper;

    private double sd = 1.0; //todo: put magic number here.
//    private ParametricMultivariateDistributionModel distribution;
    private JointBayesianBridgeDistributionModel distribution;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution acbr, double epsilon, double targetProb) {
        this.acbr = acbr;
        this.distribution = (distribution instanceof JointBayesianBridgeDistributionModel) ?
                (JointBayesianBridgeDistributionModel) distribution : null;
        this.branchRateModel = acbr.getBranchRateModel();

        if (epsilon != 0.0) {
            this.epsilon = epsilon;
        } else if (targetProb != 0.0) {
            this.targetProb = targetProb;
            this.epsilon = findEpsilonFromTargetProb(targetProb);
        }

        this.dim = acbr.getDimension();
        this.classified = new double[dim];
        classify();

        helper = new Helper();
        setupTraits();
    }

    private double findEpsilonFromTargetProb(double targetProb) {
        // finds eps such that Prob{abs(increment) < eps} = targetProb
        double probToFindEpsilon = (targetProb / 2) + 0.5;
        NormalDistribution normal = new NormalDistribution(0, sd);
        return normal.quantile(probToFindEpsilon);
    }

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

            public Double getTrait(Tree tree, NodeRef node)
            {
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

            double epsilon = xo.getAttribute(EPSILON, 0.0);
            double targetProbability = xo.getAttribute(TARGET_PROBABILITY, 0.0);

            if (epsilon < 0) {
                throw new XMLParseException("epsilon must be positive.");
            }
            else if (targetProbability < 0) {
                throw new XMLParseException("target probability must be positive.");
            }

            IncrementClassifier incrementClassifier = new IncrementClassifier(acbr, epsilon, targetProbability);
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
                AttributeRule.newDoubleRule(EPSILON, true),
                AttributeRule.newDoubleRule(TARGET_PROBABILITY, true),
        };
    };
}
