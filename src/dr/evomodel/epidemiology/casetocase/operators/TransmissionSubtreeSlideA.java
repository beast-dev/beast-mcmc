package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the subtree slide move where the slide does not affect the transmission tree.
 *
 * @author Matthew Hall
 *
 */
public class TransmissionSubtreeSlideA extends AbstractTreeOperator implements CoercableMCMCOperator {

    private CaseToCaseTreeLikelihood c2cLikelihood;
    private TreeModel tree;
    private AbstractCase[] branchMap;
    private double size = 1.0;
    private boolean gaussian = false;
    private final boolean swapInRandomRate;
    private final boolean swapInRandomTrait;
    private CoercionMode mode = CoercionMode.DEFAULT;
    private boolean debug = true;
    public static final String TRANSMISSION_SUBTREE_SLIDE_A = "transmissionSubtreeSlideA";
    public static final String SWAP_RATES = "swapInRandomRate";
    public static final String SWAP_TRAITS = "swapInRandomTrait";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";


    public TransmissionSubtreeSlideA(CaseToCaseTreeLikelihood c2cLikelihood, double weight, double size,
                                     boolean gaussian, boolean swapRates, boolean swapTraits,  CoercionMode mode) {
        this.c2cLikelihood = c2cLikelihood;
        tree = c2cLikelihood.getTree();
        branchMap = c2cLikelihood.getBranchMap();
        setWeight(weight);

        if (size == 0.0) {
            double b = 0.0;
            for (int k = 0; k < tree.getNodeCount(); ++k) {
                b += tree.getBranchLength(tree.getNode(k));
            }
            size = b / (2 * tree.getNodeCount());
        }

        this.size = size;
        this.gaussian = gaussian;
        this.swapInRandomRate = swapRates;
        this.swapInRandomTrait = swapTraits;

        this.mode = mode;
    }

    /**
     * Do a probablistic subtree slide move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() throws OperatorFailedException {


        double logq;

        final NodeRef root = tree.getRoot();
        final double oldTreeHeight = tree.getNodeHeight(root);

        NodeRef i;

        // 1. choose a random eligible node avoiding root
        do {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
        } while (root == i && !eligibleForMove(i, tree, branchMap));

        final NodeRef iP = tree.getParent(i);
        final NodeRef CiP = getOtherChild(tree, iP, i);
        final NodeRef PiP = tree.getParent(iP);

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = tree.getNodeHeight(iP);
        final double newHeight = oldHeight + delta;

        // 3. if the move is down
        if (delta > 0) {

            // 3.1 if the topology will change
            if (PiP != null && tree.getNodeHeight(PiP) < newHeight) {
                // find new parent
                NodeRef newParent = PiP;
                NodeRef newChild = iP;
                while (tree.getNodeHeight(newParent) < newHeight) {
                    newChild = newParent;
                    newParent = tree.getParent(newParent);
                    if (newParent == null) break;
                }

                // if the parent has slid out of its partition
                if(branchMap[newChild.getNumber()]!=branchMap[iP.getNumber()]){
                    return Double.NEGATIVE_INFINITY;
                }

                tree.beginTreeEdit();

                // 3.1.1 if creating a new root
                if (tree.isRoot(newChild)) {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.setRoot(iP);
                    //System.err.println("Creating new root!");

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(newChild, iP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(newChild);
                        tree.setNodeRate(newChild, tree.getNodeRate(iP));
                        tree.setNodeRate(iP, rootNodeRate);
                    }

                    // **********************************************

                }
                // 3.1.2 no new root
                else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("No new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                tree.endTreeEdit();

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(tree, newChild, oldHeight, branchMap[iP.getNumber()],
                        null);
                //System.out.println("possible sources = " + possibleSources);

                logq = -Math.log(possibleSources);

            } else {
                // just change the node height
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree up.
        else {

            // 4.0 is it a valid move?
            if (tree.getNodeHeight(i) > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (tree.getNodeHeight(CiP) > newHeight) {

                List<NodeRef> newChildren = new ArrayList<NodeRef>();
                final int possibleDestinations = intersectingEdges(tree, CiP, newHeight, branchMap[iP.getNumber()],
                        newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = MathUtils.nextInt(newChildren.size());
                NodeRef newChild = newChildren.get(childIndex);
                NodeRef newParent = tree.getParent(newChild);

                tree.beginTreeEdit();

                // 4.1.1 if iP was root
                if (tree.isRoot(iP)) {
                    // new root is CiP
                    tree.removeChild(iP, CiP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(newParent, iP);
                    tree.setRoot(CiP);

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates, so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(iP, CiP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(iP);
                        tree.setNodeRate(iP, tree.getNodeRate(CiP));
                        tree.setNodeRate(CiP, rootNodeRate);
                    }

                    // **********************************************

                    //System.err.println("DOWN: Creating new root!");
                } else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("DOWN: no new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                tree.endTreeEdit();

                logq = Math.log(possibleDestinations);
            } else {
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }

        if (swapInRandomRate) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                final double tmp = tree.getNodeRate(i);
                tree.setNodeRate(i, tree.getNodeRate(j));
                tree.setNodeRate(j, tmp);
            }

        }

        if (swapInRandomTrait) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {

                tree.swapAllTraits(i, j);

//                final double tmp = tree.getNodeTrait(i, TRAIT);
//                tree.setNodeTrait(i, TRAIT, tree.getNodeTrait(j, TRAIT));
//                tree.setNodeTrait(j, TRAIT, tmp);
            }

        }

        if (logq == Double.NEGATIVE_INFINITY) throw new OperatorFailedException("invalid slide");

        if (debug) c2cLikelihood.checkPartitions();

        return logq;
    }

    private double getDelta() {
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
    }

    private boolean eligibleForMove(NodeRef node, TreeModel tree, AbstractCase[] branchMap){
        // to be eligible for this move, the node's parent and grandparent, or parent and other child, must be in the
        // same partition (so removing the parent has no effect on the transmission tree)
        // if only one of these things is true you can only slide one way - might be able to do better than this

        return branchMap[tree.getParent(node).getNumber()]==branchMap[tree.getParent(tree.getParent(node)).getNumber()]
                || branchMap[tree.getParent(node).getNumber()]==branchMap[getOtherChild(tree,
                tree.getParent(node), node).getNumber()];
    }

    //intersectingEdges here is modified to count only possible sources for this special case of the operator - i.e.
    //only branches which have one end in the relevant partition

    private int intersectingEdges(Tree tree, NodeRef node, double height, AbstractCase partition,
                                  List<NodeRef> directChildren) {

        final NodeRef parent = tree.getParent(node);

        if (tree.getNodeHeight(parent) < height || branchMap[parent.getNumber()]!=partition) return 0;

        if (tree.getNodeHeight(node) < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += intersectingEdges(tree, tree.getChild(node, i), height, partition, directChildren);
        }
        return count;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    public double getCoercableParameter() {
        return Math.log(getSize());
    }

    public void setCoercableParameter(double value) {
        setSize(Math.exp(value));
    }

    public double getRawParameter() {
        return getSize();
    }

    public CoercionMode getMode() {
        return mode;
    }


    public String getPerformanceSuggestion() {
        return "not implemented";
    }

    public String getOperatorName() {
        return TRANSMISSION_SUBTREE_SLIDE_A + "(" + tree.getId() + ")";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean swapRates = xo.getAttribute(SWAP_RATES, false);
        boolean swapTraits = xo.getAttribute(SWAP_TRAITS, false);

        CoercionMode mode = CoercionMode.DEFAULT;
        if (xo.hasAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
            if (xo.getBooleanAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                mode = CoercionMode.COERCION_ON;
            } else {
                mode = CoercionMode.COERCION_OFF;
            }
        }

        CaseToCaseTreeLikelihood c2cL = (CaseToCaseTreeLikelihood)xo.getChild(CaseToCaseTreeLikelihood.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        final double size = xo.getAttribute("size", 1.0);

        if (Double.isInfinite(size) || size <= 0.0) {
            throw new XMLParseException("size attribute must be positive and not infinite. was " + size +
                    " for tree " + c2cL.getTree().getId() );
        }

        final boolean gaussian = xo.getBooleanAttribute("gaussian");
        TransmissionSubtreeSlideA operator = new TransmissionSubtreeSlideA(c2cL, weight, size, gaussian,
                swapRates, swapTraits, mode);
        operator.setTargetAcceptanceProbability(targetAcceptance);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that slides a phylogenetic subtree, preserving the transmission tree topology.";
    }

    public Class getReturnType() {
        return TransmissionSubtreeSlideA.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            // Make size optional. If not given or equals zero, size is set to half of average tree branch length.
            AttributeRule.newDoubleRule("size", true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule("gaussian"),
            AttributeRule.newBooleanRule(SWAP_RATES, true),
            AttributeRule.newBooleanRule(SWAP_TRAITS, true),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}
