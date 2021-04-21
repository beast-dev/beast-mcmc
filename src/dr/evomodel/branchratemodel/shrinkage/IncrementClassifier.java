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
    public static final String BY_SIGN = "bySign";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";


    private AutoCorrelatedBranchRatesDistribution acbr; //autocorrelated branch rates
    private DifferentiableBranchRates branchRateModel;
    private final classificationMode classifier;
    private double epsilon;
    private int dim;
    private Helper helper;

    private double[] classified;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution acbr, double epsilon, boolean bySign) {
        this.acbr = acbr;
        this.branchRateModel = acbr.getBranchRateModel();
        this.epsilon = epsilon;
        this.classifier = bySign ? classificationMode.BY_SIGN : classificationMode.BY_EPSILON;

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
            classified[i] = 0;
            if (classifier.getIncrement(increment) > epsilon) {
                classified[i] = 1;
            }
        }
    }

    public enum classificationMode {

        BY_EPSILON("byEpsilon") {
            @Override
            double getIncrement(double increment) {
                return (Math.abs(increment));
            }
        },

        BY_SIGN("bySign") {
            @Override
            double getIncrement(double increment) {
                return increment;
            }
        };

        classificationMode(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        private final String name;

        abstract double getIncrement(double increment);
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

            boolean bySign = xo.getAttribute(BY_SIGN, false);

            if (bySign) {
                if (epsilon != 0.0) {
                    throw new RuntimeException("Sign classifier should use epsilon = 0.0.");
                }
            } else if (!bySign) {
                if (!xo.hasAttribute(EPSILON) && !xo.hasAttribute(TARGET_PROBABILITY)) {
                    throw new RuntimeException("Must specify epsilon or target probability when not using sign classifier.");
                }
            }

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

            IncrementClassifier incrementClassifier = new IncrementClassifier(acbr, epsilon, bySign);
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
            return "Classifies increment as 0 or 1 based on sign or, alternatively, arbitrary cutoff epsilon.";
        }

        @Override
        public Class getReturnType() {
            return IncrementClassifier.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(AutoCorrelatedBranchRatesDistribution.class),
//                AttributeRule.newStringRule(CLASSIFICATION_MODE, false),
                new XORRule(
                        AttributeRule.newBooleanRule(BY_SIGN, false),
                        new XORRule(
                                AttributeRule.newDoubleRule(EPSILON, false),
                                new AndRule(
                                        AttributeRule.newDoubleRule(TARGET_PROBABILITY, false),
                                        new ElementRule(BayesianBridgePriorSampler.class, false)
                                )
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