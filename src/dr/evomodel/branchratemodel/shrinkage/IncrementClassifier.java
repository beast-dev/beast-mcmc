package dr.evomodel.branchratemodel.shrinkage;


import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class IncrementClassifier implements TreeTraitProvider {

    public static final String PARAMETER_CLASSIFIER = "parameterClassifier";
    public static final String EPSILON = "epsilon";
    public static final String TARGET_PROBABILITY = "targetProbability";

    private AutoCorrelatedBranchRatesDistribution brm;
    private double epsilon;
    private double targetProb;

    public IncrementClassifier(AutoCorrelatedBranchRatesDistribution brm, double epsilon, double targetProb) {
        this.brm = brm;
        if (epsilon != 0.0){
            this.epsilon = epsilon;
        }
        else if (targetProb != 0.0){
            this.targetProb = targetProb;
            this.epsilon = findEpsilonFromTargetProb(targetProb);
        }
    }

    private double findEpsilonFromTargetProb(double targetProb){
        return 0.0;
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

            if (xo.hasAttribute(EPSILON) && xo.hasAttribute(TARGET_PROBABILITY)){
                throw new XMLParseException("Cannot set both epsilon and target probability.");
            }

            double epsilon = xo.getAttribute(EPSILON, 0);
            double targetProbability = xo.getAttribute(TARGET_PROBABILITY, 0);

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
