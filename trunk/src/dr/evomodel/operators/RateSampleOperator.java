package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.GibbsOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * A special operator for sampling rates in a subtree
 * according to a RateEvolutionLikelihood such as the autocorrelated clock model
 *
 * @author Michael Defoin Platel
 */
public class RateSampleOperator extends SimpleMCMCOperator {

    public static final String SAMPLE_OPERATOR = "rateSampleOperator";
    public static final String SAMPLE_ALL = "sampleAll";

    private TreeModel tree;

    private boolean sampleAll;

    RateEvolutionLikelihood rateEvolution;

    public RateSampleOperator(TreeModel tree, boolean sampleAll, RateEvolutionLikelihood rateEvolution) {

        this.tree = tree;
        this.sampleAll = sampleAll;
        this.rateEvolution = rateEvolution;
    }

    /**
     * sample the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        int index;
        if (sampleAll) {
            index = tree.getRoot().getNumber();
        } else {
            do {
                index = MathUtils.nextInt(tree.getNodeCount());
            } while (tree.isExternal(tree.getNode(index)));
        }


        double logBackward = rateEvolution.getLogLikelihood();

        //sampleOne(tree.getNode(index));

        //sampleSubtree(tree.getNode(index));
        sampleNode(tree.getNode(index));

        //sampleSister(tree.getNode(index));

        double logForward = rateEvolution.getLogLikelihood();

        return logBackward - logForward;
    }

    void sampleSubtree(NodeRef parent) {

        int nbChildren = tree.getChildCount(parent);
        for (int c = 0; c < nbChildren; c++) {
            final NodeRef node = tree.getChild(parent, c);
            rateEvolution.sampleRate(node);
            sampleSubtree(node);
        }
    }

    void sampleSister(NodeRef parent) {

        int nbChildren = tree.getChildCount(parent);
        for (int c = 0; c < nbChildren; c++) {
            final NodeRef node = tree.getChild(parent, c);
            rateEvolution.sampleRate(node);
        }

    }

    void sampleNode(NodeRef parent) {

        int nbChildren = tree.getChildCount(parent);

        if (nbChildren > 0) {
            final int c = MathUtils.nextInt(nbChildren);

            final NodeRef node = tree.getChild(parent, c);
            rateEvolution.sampleRate(node);
        }

    }

    void sampleOne(NodeRef parent) {

        int nbChildren = tree.getChildCount(parent);

        if (nbChildren > 0) {

            final NodeRef node = tree.getChild(parent, 0);
            rateEvolution.sampleRate(node);
        }

    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "rateSample";
    }


    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        return "No suggestions";
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SAMPLE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            final double weight = xo.getDoubleAttribute(WEIGHT);

            final boolean sampleAll = xo.getBooleanAttribute(SAMPLE_ALL);

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            RateEvolutionLikelihood rateEvolution = (RateEvolutionLikelihood) xo.getChild(RateEvolutionLikelihood.class);

            RateSampleOperator operator = new RateSampleOperator(treeModel, sampleAll, rateEvolution);
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a rateSample operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(SAMPLE_ALL, true),
                new ElementRule(TreeModel.class),
                new ElementRule(RateEvolutionLikelihood.class, true),
        };

    };

    public String toString() {
        return "rateSampleOperator(";
    }

    //PRIVATE STUFF

}