package dr.evomodel.operators;

import dr.inference.operators.*;
import dr.inference.model.Parameter;
import dr.inference.model.Bounds;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import dr.evolution.tree.NodeRef;
import dr.xml.*;

import java.util.List;
import java.util.ArrayList;

/**
 * A special operator for scaling the variance of the autocorrelated clock model
 * and subsequnetly the rates in a tree.
 *
 * @author Michael Defoin Platel
 */
public class RateVarianceScaleOperator extends AbstractCoercableOperator {

    public static final String SCALE_OPERATOR = "rateVarianceScaleOperator";
    public static final String SCALE_FACTOR = "scaleFactor";

    private TreeModel tree;
    private Parameter variance;

    public RateVarianceScaleOperator(TreeModel tree, Parameter variance, double scale, CoercionMode mode) {
        super(mode);

        this.scaleFactor = scale;
        this.tree = tree;
        this.variance = variance;
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        //Scale the variance
        double oldValue = variance.getParameterValue(0);
        double newValue = scale * oldValue;
        double logq= - Math.log(scale);

        final Bounds bounds = variance.getBounds();
        if (newValue < bounds.getLowerLimit(0) || newValue > bounds.getUpperLimit(0)) {
            throw new OperatorFailedException("proposed value outside boundaries");
        }
        variance.setParameterValue(0, newValue);

        //Scale the rates of the tree accordingly
        NodeRef root =  tree.getRoot();
        final int index = root.getNumber();

        List<NodeRef> listNode = new ArrayList<NodeRef>();
        getSubtree(listNode, tree.getNode(index));

        final double rateScale = Math.sqrt(scale);

        for( NodeRef node : listNode){

            oldValue = tree.getNodeRate(node);
            newValue = oldValue * rateScale;

            tree.setNodeRate(node, newValue);
        }

        //  According to the hastings ratio in the scale Operator
        logq += (listNode.size()  - 2) * Math.log(rateScale);

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
        return "rateVarianceScale";
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

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            final Parameter variance = (Parameter) xo.getChild(Parameter.class);
            if (variance.getDimension() != 1) {
                throw new XMLParseException("dimension of the variance parameter should be 1");
            }

            RateVarianceScaleOperator operator = new RateVarianceScaleOperator(treeModel, variance, scaleFactor, mode);
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
                new ElementRule(Parameter.class),
        };

    };

    public String toString() {
        return "rateVarianceScaleOperator(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}
