package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Alexander Fisher
 * Symmetric rooted NNI for use with marginal tree likelihood via ApproximateTreeDataLikelihood
 *
 * For n tips, choose one of n-2 non-root internal nodes uniformly, then one of
 * its two children uniformly. Swap that child with the node's sibling. Each
 * topology neighbor has probability 1/[2(n-2)], including its reverse move.
 * Neither selection nor feasibility depends on internal heights.
 *
 * ONLY for a collapsed topology target: internal heights are computational
 * workspace, not sampled state variables. After each swap they are reset
 * deterministically, in postorder, to max(child heights) + heightStep. Tip
 * dates are preserved. The likelihood must integrate out internal heights
 * (e.g. Laplace evidence) and produce a topology-dependent score independent
 * of incoming workspace heights. There is no continuous-height Hastings term.
 * This operator is NOT valid with an ordinary joint tree/height posterior.
 *
 * No monophyly, calibration, or fixed-height constraints are enforced here.
 * Hard topology constraints should reject the proposed state in the target;
 * do not resample until an allowed topology is obtained.
 */
public final class TopologyOnlyNNI extends AbstractTreeOperator {
    public static final String NAME = "topologyOnlyNNI";
    private final TreeModel tree;
    private final double heightStep;

    public TopologyOnlyNNI(TreeModel tree, double weight) {
        this(tree, weight, 1.0);
    }

    public TopologyOnlyNNI(TreeModel tree, double weight, double heightStep) {
        if (tree instanceof ConstrainedTreeModel) {
            throw new IllegalArgumentException("topologyOnlyNNI requires an unconstrained tree model");
        }
        if (!(weight > 0.0) || !Double.isFinite(weight)) {
            throw new IllegalArgumentException("weight must be finite and positive");
        }
        if (!(heightStep > 0.0) || !Double.isFinite(heightStep)) {
            throw new IllegalArgumentException("heightStep must be finite and positive");
        }
        if (tree.getExternalNodeCount() < 3) {
            throw new IllegalArgumentException("topologyOnlyNNI needs at least three tips");
        }
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            if (tree.getChildCount(tree.getInternalNode(i)) != 2) {
                throw new IllegalArgumentException("topologyOnlyNNI requires a strictly binary tree");
            }
        }
        this.tree = tree;
        this.heightStep = heightStep;
        setWeight(weight);
        validateHeightRange();
    }

    @Override
    public double doOperation() {
        return propose(MathUtils.nextInt(tree.getInternalNodeCount() - 1), MathUtils.nextInt(2));
    }

    // Package-private deterministic choice for exhaustive transition tests.
    double propose(int edgeIndex, int childIndex) {
        if (edgeIndex < 0 || edgeIndex >= tree.getInternalNodeCount() - 1
                || childIndex < 0 || childIndex > 1) {
            throw new IllegalArgumentException("Invalid NNI choice");
        }
        validateHeightRange(); // check before mutating, including changed tip dates
        NodeRef childParent = null;
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef candidate = tree.getInternalNode(i);
            if (!tree.isRoot(candidate) && edgeIndex-- == 0) {
                childParent = candidate;
                break;
            }
        }
        NodeRef parent = tree.getParent(childParent);
        NodeRef child = tree.getChild(childParent, childIndex);
        NodeRef uncle = getOtherChild(tree, parent, childParent);

        tree.beginTreeEdit();
        try {
            tree.removeChild(childParent, child);
            tree.removeChild(parent, uncle);
            tree.addChild(childParent, uncle);
            tree.addChild(parent, child);
            // Reset before endTreeEdit dispatches events: observers see a valid tree.
            resetHeights(tree.getRoot());
        } finally {
            tree.endTreeEdit();
        }
        return 0.0;
    }

    private double resetHeights(NodeRef node) {
        if (tree.isExternal(node)) {
            return tree.getNodeHeight(node);
        }
        double lower = Math.max(resetHeights(tree.getChild(node, 0)),
                                resetHeights(tree.getChild(node, 1)));
        double height = nextHeight(lower);
        tree.setNodeHeight(node, height);
        return height;
    }

    private double nextHeight(double lower) {
        double height = lower + heightStep;
        return height > lower ? height : Math.nextUp(lower);
    }

    private void validateHeightRange() {
        double upper = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            double height = tree.getNodeHeight(tree.getExternalNode(i));
            if (!Double.isFinite(height)) {
                throw new IllegalArgumentException("Tip heights must be finite");
            }
            upper = Math.max(upper, height);
        }
        // A caterpillar is the deepest possible binary topology.
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            upper = nextHeight(upper);
        }
        if (!Double.isFinite(upper)) {
            throw new IllegalArgumentException("Tip heights/heightStep overflow internal heights");
        }
    }

    @Override
    public String getOperatorName() {
        return NAME + "(" + tree.getId() + ")";
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public String getParserName() { return NAME; }
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            try {
                return new TopologyOnlyNNI((TreeModel) xo.getChild(TreeModel.class),
                        xo.getDoubleAttribute(MCMCOperator.WEIGHT), xo.getAttribute("heightStep", 1.0));
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(e.getMessage());
            }
        }
        public String getParserDescription() {
            return "Symmetric topology-only rooted NNI for a collapsed target; resets internal heights.";
        }
        public Class getReturnType() { return TopologyOnlyNNI.class; }
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule("heightStep", true),
                new ElementRule(TreeModel.class)
            };
        }
    };
}
