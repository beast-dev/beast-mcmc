package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator swaps a randomly selected rate change from parent to offspring or vice versa.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeBitMoveOperator extends SimpleMCMCOperator {

    public static final String BIT_MOVE_OPERATOR = "treeBitMoveOperator";
    public static final String INDICTATOR_TRAIT = "indicatorTrait";
    public static final String TRAIT2 = "trait2";


    public TreeBitMoveOperator(TreeModel tree, String t1, String t2, double weight) {
        this.tree = tree;
        this.indicatorTrait = t1;
        this.trait2 = t2;

        if (indicatorTrait == null) indicatorTrait = "trait";

        setWeight(weight);
    }

    /**
     * Pick a parent-child node pair involving a single rate change and swap the rate change location
     * and corresponding rate parameters.
     */
    public final double doOperation() throws OperatorFailedException {

        NodeRef root = tree.getRoot();

        // 1. collect nodes that form a pair with parent such that
        // one of them has a one and one has a zero
        List<NodeRef> candidates = new ArrayList<NodeRef>();

        for (int i = 0; i < tree.getNodeCount(); i++) {

            NodeRef node = tree.getNode(i);
            if (node != root && tree.getParent(node) != root) {

                NodeRef parent = tree.getParent(node);

                int sum = rateChange(tree, node) + rateChange(tree, parent);

                if (sum == 1) candidates.add(node);
            }
        }

        if (candidates.size() == 0) throw new OperatorFailedException("No suitable pairs!");

        NodeRef node = candidates.get(MathUtils.nextInt(candidates.size()));
        NodeRef parent = tree.getParent(node);

        double nodeTrait, parentTrait;
        double nodeRate, parentRate;

        nodeTrait = tree.getNodeTrait(node, indicatorTrait);
        parentTrait = tree.getNodeTrait(parent, indicatorTrait);

        tree.setNodeTrait(node, indicatorTrait, parentTrait);
        tree.setNodeTrait(parent, indicatorTrait, nodeTrait);

        if (trait2 != null) {
            nodeTrait = tree.getNodeTrait(node, trait2);
            parentTrait = tree.getNodeTrait(parent, trait2);

            tree.setNodeTrait(node, trait2, parentTrait);
            tree.setNodeTrait(parent, trait2, nodeTrait);
        } else {
            nodeRate = tree.getNodeRate(node);
            parentRate = tree.getNodeRate(parent);

            tree.setNodeRate(node, parentRate);
            tree.setNodeRate(parent, nodeRate);
        }

        return 0.0;
    }

    public final int rateChange(TreeModel tree, NodeRef node) {
        return (int) Math.round(tree.getNodeTrait(node, indicatorTrait));
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "treeBitMove()";
    }

    public final String getPerformanceSuggestion() {

        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIT_MOVE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);


            String trait1 = null;
            String trait2 = null;
            if (xo.hasAttribute(INDICTATOR_TRAIT)) trait1 = xo.getStringAttribute(INDICTATOR_TRAIT);
            if (xo.hasAttribute(TRAIT2)) trait2 = xo.getStringAttribute(TRAIT2);

            return new TreeBitMoveOperator(treeModel, trait1, trait2, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a bit-move operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(TreeModel.class),
                AttributeRule.newStringRule(INDICTATOR_TRAIT, true),
                AttributeRule.newStringRule(TRAIT2, true)
        };

    };
    // Private instance variables

    private TreeModel tree;
    private String indicatorTrait;
    private String trait2;
}
