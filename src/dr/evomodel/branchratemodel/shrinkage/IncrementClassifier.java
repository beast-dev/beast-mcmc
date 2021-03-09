package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.operators.shrinkage.BayesianBridgePriorSampler;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider, Loggable {
    public static final String INCREMENT_CLASSIFIER = "incrementClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";

    private AutoCorrelatedBranchRatesDistribution acbr; //autocorrelated branch rates
    private DifferentiableBranchRates branchRateModel;
    private double epsilon;
    private int dim;
    private Helper helper;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution acbr, double epsilon) {
        this.acbr = acbr;
        this.branchRateModel = acbr.getBranchRateModel();
        this.epsilon = epsilon;
        this.dim = acbr.getDimension();
        this.classified = new double[dim];
        classify();
        System.out.println("EPSILON is " + epsilon);
        helper = new Helper();
        setupTraits();
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

            AutoCorrelatedBranchRatesDistribution acbr = (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);

            double epsilon = xo.getAttribute(EPSILON, 0.0);
            double targetProbability = xo.getAttribute(TARGET_PROBABILITY, 0.0);

            if (epsilon < 0.0) {
                throw new XMLParseException("epsilon must be positive.");
            } else if (targetProbability < 0.0) {
                throw new XMLParseException("target probability must be positive.");
            }

            BayesianBridgePriorSampler bayesianBridgePriorSampler;
            if (targetProbability != 0.0) {
                bayesianBridgePriorSampler = (BayesianBridgePriorSampler) xo.getChild(BayesianBridgePriorSampler.class);
                double steps = bayesianBridgePriorSampler.getSteps();
                // ensure at least 50 samples from median, in general: 2*n / steps for n steps away
                if (targetProbability < 100 / steps) {
                    throw new XMLParseException("For BB prior with " + bayesianBridgePriorSampler.getSteps() + " steps, target Probability should be greater than " + 100.0 / steps);
                }
                epsilon = bayesianBridgePriorSampler.getEpsilon(targetProbability);
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
        LogColumn[] columns = new LogColumn[dim];
        classify();
        for (int i = 0; i < dim; i++) {
            String colName = "incrementClass.";
            int finalI = i;
            columns[i] = new NumberColumn(colName + finalI) {
                @Override
                public double getDoubleValue() {
                        return classified[finalI];
                }
            };
        }
        return columns;
    }
}