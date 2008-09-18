package dr.evomodel.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.NodeRef;

import java.util.List;
import java.util.ArrayList;

/**
 * A special operator for scaling rates in a subtree.
 *
 * @author Michael Defoin Platel
 */
public class RateScaleOperator extends AbstractCoercableOperator {

    public static final String SCALE_OPERATOR = "rateScaleOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String NO_ROOT = "noRoot";

    private TreeModel tree;

    private boolean noRoot;

    public RateScaleOperator(TreeModel tree, double scale, boolean noRoot, CoercionMode mode) {

        super(mode);

        this.scaleFactor = scale;
        this.tree = tree;

        this.noRoot = noRoot;
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        int index;
        if (noRoot){
        do{
            index = MathUtils.nextInt(tree.getNodeCount());
        }while(index == tree.getRoot().getNumber());
        }else{
            index = MathUtils.nextInt(tree.getNodeCount());
        }

        //NodeRef root =  tree.getRoot();
        //index = root.getNumber(); 

        List<NodeRef> listNode = new ArrayList<NodeRef>();
        getSubtree(listNode, tree.getNode(index));

        double oldValue, newValue;
        double logq=0;
        for( NodeRef node : listNode){

            oldValue = tree.getNodeRate(node);
            newValue = oldValue * scale;

            tree.setNodeRate(node, newValue);
        }

        //  According to the hastings ratio in the scale Operator
        logq = (listNode.size()  - 2) * Math.log(scale);

        return logq;
    }

    void getSubtree(List<NodeRef> listNode, NodeRef parent){

        listNode.add(parent);
        int nbChildren = tree.getChildCount(parent);
        for (int c=0; c < nbChildren; c++ ){
            getSubtree(listNode, tree.getChild(parent, c));
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
        return "rateScale";
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
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

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return SCALE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

           CoercionMode mode = CoercionMode.parseMode(xo);

            final double weight = xo.getDoubleAttribute(WEIGHT);
            final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            final boolean noRoot = xo.getBooleanAttribute(NO_ROOT);

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            RateScaleOperator operator = new RateScaleOperator(treeModel, scaleFactor, noRoot, mode);
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a rateScale operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
				new ElementRule(TreeModel.class),
        };

    };

    public String toString() {
        return "rateScaleOperator(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}
