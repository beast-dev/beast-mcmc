/*
 * ObsoleteARGNewEventOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/*
 * AddRemoveSubtreeOperator.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.evomodel.arg.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGModel.Node;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.functionEval.GammaFunction;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Implements the subtree slide move.
 *
 * @author Marc Suchard
 * @version $Id: ObsoleteARGAddRemoveEventOperator.java,v 1.18.2.4 2006/11/06 01:38:30 msuchard Exp $
 */

public class ObsoleteARGNewEventOperator extends AbstractCoercableOperator {
//		SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String SUBTREE_SLIDE = "newARGEvent";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_TRAITS = "swapTraits";
    public static final String MAX_VALUE = "maxTips";
    public static final String SINGLE_PARTITION = "singlePartitionProbability";
    public static final String IS_RECOMBINATION = "isRecombination";

    public static final String JUST_INTERNAL = "justInternalNodes";
    public static final String INTERNAL_AND_ROOT = "internalAndRootNodes";
    public static final String NODE_RATES = TreeModelParser.NODE_RATES;

    private ARGModel arg = null;
    private double size = 1.0;
    private boolean gaussian = false;
    private double singlePartitionProbability = 0.0;
    private boolean isRecombination = false;
    //   private boolean swapRates;
    //   private boolean swapTraits;
    //	private int mode = CoercableMCMCOperator.DEFAULT;
    private CompoundParameter internalNodeParameters;
    private CompoundParameter internalAndRootNodeParameters;
    private CompoundParameter nodeRates;
//	private int maxTips = 1;

    public ObsoleteARGNewEventOperator(ARGModel arg, int weight, double size, boolean gaussian,
                                       boolean swapRates, boolean swapTraits, CoercionMode mode,
                                       CompoundParameter param1,
                                       CompoundParameter param2,
                                       CompoundParameter param3,
                                       double singlePartitionProbability, boolean isRecombination) {
        super(mode);
        this.arg = arg;
        setWeight(weight);
//		this.maxTips = maxTips;
        this.size = size;
        this.gaussian = gaussian;
//        this.swapRates = swapRates;
//        this.swapTraits = swapTraits;
        this.internalNodeParameters = param1;
        this.internalAndRootNodeParameters = param2;
        this.nodeRates = param3;
        this.singlePartitionProbability = singlePartitionProbability;
        this.isRecombination = isRecombination;

//		this.mode = mode;
    }


    /**
     * Do a add/remove reassortment node operation
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
//		System.err.println("Starting AddRemove Operation");

        double logq = 0;
        try {
            if (MathUtils.nextDouble() < 0.5)
                logq = AddOperation();
            else
                logq = RemoveOperation();
        } catch (Exception ofe) {
//			System.err.println(ofe);
            if (ofe.getMessage().compareTo("No reassortment nodes to remove.") != 0) {
                System.err.println("Catch: " + ofe.getMessage());
                System.exit(-1);
            }
        }
//        if (arg.isBifurcationDoublyLinked(arg.getRoot()))
//            throw new OperatorFailedException("trouble with double-rooted root");
//	    System.err.println("logq = "+logq);
        return logq;
    }

    private int findPotentialNodesToRemove(ArrayList<NodeRef> list) {
        int count = 0;
        int n = arg.getNodeCount();
//	    int max = arg.getReassortmentNodeCount();
        Node root = (Node) arg.getRoot();
        for (int i = 0; i < n; i++) {
            Node node = (Node) arg.getNode(i);
            if (node.isReassortment()
                    && ((node.leftParent != root && node.leftParent.isBifurcation()) ||
                    (node.rightParent != root && node.rightParent.isBifurcation()))
//	                && !(
//		                    (node.leftParent == root && node.rightParent.isReassortment()) ||
//				            (node.rightParent== root && node.leftParent.isReassortment()))
//		            && !(node.leftParent == root && node.rightParent == root)
//		            && (node.leftParent.isBifurcation() || node.rightParent.isBifurcation())
                    ) {
                if (list != null)
                    list.add(node);
                count++;
            }
        }
//	    System.err.printf("max = %d, available = %d\n",max,count);
        return count;
    }

    private double RemoveOperation() throws ARGOperatorFailedException {
        double logq = 0;

//	    System.err.println("Starting remove ARG operation.");
//	    System.err.println("Before remove:\n"+arg.toGraphString());

        // 1. Draw reassortment node uniform randomly

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();

        int totalPotentials = findPotentialNodesToRemove(potentialNodes);
        if (totalPotentials == 0)
            throw new ARGOperatorFailedException("No reassortment nodes to remove.");
//	    System.err.println("potentials exist!");
        Node recNode = (Node) potentialNodes.get(MathUtils.nextInt(totalPotentials));
//        logq += Math.log(totalPotentials);

        double reverseReassortmentHeight = 0;
        double reverseBifurcationHeight = 0;

        double treeHeight = arg.getNodeHeight(arg.getRoot());

//	    System.err.println("RecNode to remove = "+recNode.number);
        reverseReassortmentHeight = arg.getNodeHeight(recNode);

        arg.beginTreeEdit();
        boolean doneSomething = false;
        Node recParent = recNode.leftParent;
        Node recChild = recNode.leftChild;
        if (recNode.leftParent == recNode.rightParent) { // Doubly linked.
            Node recGrandParent = recParent.leftParent;
//	        reverseBifurcationHeight = arg.getNodeHeight(recParent);

//            reverseReassortmentSpan = arg.getNodeHeight(recParent) - arg.getNodeHeight(recChild);
//            reverseBifurcationSpan = arg.getNodeHeight(recGrandParent) - arg.getNodeHeight(recChild);
//            if (arg.isRoot(recParent)) { // This case should never happen as double links
//                arg.setRoot(recChild);    // to root can not be added or removed.
//	            System.err.println("doubled link root?");
//	            System.exit(-1);
//            } else {
            //Node recGrandParent = recParent1.leftParent; // And currently recParent must be a bifurcatio
            //System.err.println("recGrand   ="+recGrandParent.number);
            arg.doubleRemoveChild(recGrandParent, recParent);
            arg.doubleRemoveChild(recNode, recChild);
            if (recGrandParent.bifurcation)
                arg.singleAddChild(recGrandParent, recChild);
            else
                arg.doubleAddChild(recGrandParent, recChild);
//            }
            doneSomething = true;
            // There are not left/right choices to be made for doubly linked removals.
            // End doubly linked.
        } else { // Two different parents.
            Node recParent1 = recNode.leftParent;
            Node recParent2 = recNode.rightParent;
            if (!recParent1.bifurcation || recParent1.isRoot()) { // One orientation is valid
                recParent1 = recNode.rightParent;
                recParent2 = recNode.leftParent;
            } else if (recParent2.bifurcation && !recParent2.isRoot()) { // Both orientations are valid
                if (MathUtils.nextDouble() < 0.5) { // choose equally likely
                    recParent1 = recNode.rightParent;
                    recParent2 = recNode.leftParent;
                }
                logq += Math.log(2);
            }

            Node recGrandParent = recParent1.leftParent; // And currently recParent must be a bifurcatio

            Node otherChild = recParent1.leftChild;
            if (otherChild == recNode)
                otherChild = recParent1.rightChild;

            if (recGrandParent.bifurcation)
                arg.singleRemoveChild(recGrandParent, recParent1);
            else
                arg.doubleRemoveChild(recGrandParent, recParent1);

            arg.singleRemoveChild(recParent1, otherChild);
            if (recParent2.bifurcation)
                arg.singleRemoveChild(recParent2, recNode);
            else
                arg.doubleRemoveChild(recParent2, recNode);
            arg.doubleRemoveChild(recNode, recChild);
            if (otherChild != recChild) {
                if (recGrandParent.bifurcation)
                    arg.singleAddChild(recGrandParent, otherChild);
                else
                    arg.doubleAddChild(recGrandParent, otherChild);
                if (recParent2.bifurcation)
                    arg.singleAddChild(recParent2, recChild);
                else
                    arg.doubleAddChild(recParent2, recChild);
            } else {
                if (recGrandParent.bifurcation)
                    arg.singleAddChildWithOneParent(recGrandParent, otherChild);
                else
                    arg.doubleAddChildWithOneParent(recGrandParent, otherChild);
                if (recParent2.bifurcation)
                    arg.singleAddChildWithOneParent(recParent2, recChild);
                else
                    arg.doubleAddChildWithOneParent(recParent2, recChild);
            }
//			System.err.println("Sanity check in Remove Operator");
//			sanityCheck();
//			System.err.println("End Remove Operator in sanity check");
            doneSomething = true;

            recParent = recParent1;
        }

        reverseBifurcationHeight = arg.getNodeHeight(recParent);

        if (doneSomething) {
            try {
                arg.contractARGWithRecombinant(recParent, recNode,
                        internalNodeParameters, internalAndRootNodeParameters, nodeRates);
            } catch (Exception e) {
                System.err.println("here");
                System.err.println(e);

            }
        }

        arg.pushTreeSizeChangedEvent();
        arg.endTreeEdit();
//        try {
//            arg.checkTreeIsValid();
//        } catch (MutableTree.InvalidTreeException ite) {
//            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
//                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
//        }

//	    double d1 = findPotentialAttachmentPoints(reverseBifurcationHeight,null);
//	    double d2 = findPotentialAttachmentPoints(reverseReassortmentHeight,null);
//	    System.err.printf("d1 = %5.4f, d2 = %5.4f\n",d1,d2);

//	    logq -= Math.log(findPotentialAttachmentPoints(reverseBifurcationHeight,null));
//	    logq -= Math.log(findPotentialAttachmentPoints(reverseReassortmentHeight,null));

//        int nodes = arg.getInternalNodeCount() - 1;
//
//        logq += 3*Math.log(2);
//
//        logq = 100;

//	    System.err.println("logq remove = "+logq);
        logq = 0;
//	    logq -= 10;
//	    logq += Math.log(totalPotentials);
        return logq; // 1 / total potentials * 1 / 2 (if valid) * length1 * length2 * attachmentSisters
    }

    private static double lnGamma(double x) {
        if (x == 1 || x == 2)
            return 0.0;
        return GammaFunction.logGamma(x);
    }

    private void checkAllHeights() {
        int len = arg.getInternalNodeCount();
        System.err.println("# internal nodes = " + len);
        int n = internalNodeParameters.getParameterCount();
        System.err.println("VSCP (" + n + ")");
        for (int i = 0; i < n; i++) {
            System.err.println(internalNodeParameters.getParameterValue(i));
        }
        n = arg.getInternalNodeCount();
        System.err.println("Checking all internal nodes (" + n + ") via tree:");
        for (int i = 0; i < n; i++) {
            NodeRef node = arg.getInternalNode(i);
            System.err.print(TreeUtils.uniqueNewick(arg, node) + " ");
            System.err.println(((Node) node).getHeight());
        }
    }

    private int findPotentialAttachmentPoints(double time, ArrayList<NodeRef> list) {
        int count = 0;
        for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
            NodeRef nr = arg.getNode(i);
            if (!arg.isRoot(nr) && arg.getNodeHeight(nr) < time) {
                if (arg.getNodeHeight(arg.getParent(nr, 0)) > time) {
                    if (list != null)
                        list.add(nr);
                    count++;
                }
                if (arg.isReassortment(nr) && arg.getNodeHeight(arg.getParent(nr, 1)) > time) {
                    if (list != null)
                        list.add(nr);
                    count++;
                }
            }
        }
        return count;
    }

    private double drawRandomPartitioning(Parameter partitioning) {
        double logq = 0;
        int len = arg.getNumberOfPartitions();
        if (len == 2) {
//            boolean first = MathUtils.nextBoolean();
            if (partitioning != null) {
                if (MathUtils.nextBoolean())
                    partitioning.setParameterValueQuietly(0, 1.0);
                else
                    partitioning.setParameterValueQuietly(1, 1.0);
            }
            return Math.log(2);
        }
        if (isRecombination) {
            logq += drawRandomRecombination(partitioning);
        } else {
            logq += drawRandomReassortment(partitioning);
        }
        return logq;
    }


    /* Draws a new partitioning.
      * With probability singlePartitionProbability, one bit is set;
      * otherwise, all bits are selected via a random permutation
      *
      */
    private double drawRandomReassortment(Parameter partitioning) {
        int len = arg.getNumberOfPartitions();
        double logq = 0;
        if (MathUtils.nextDouble() < singlePartitionProbability) {
            if (partitioning != null)
                partitioning.setParameterValueQuietly(MathUtils.nextInt(len), 1.0);
            return Math.log(len);
        }
        int[] permutation = MathUtils.permuted(len);
        int cut = MathUtils.nextInt(len - 1);
        for (int i = 0; i < len; i++) {
            logq += Math.log(i + 1);
            if (i > cut && partitioning != null)
                partitioning.setParameterValueQuietly(permutation[i], 1.0);
        }
        logq += Math.log(len - 1);
        return logq;
    }

    /* Draws a new partitioning.
      * A break-pt is drawn uniformly
      *
      */
    private double drawRandomRecombination(Parameter partitioning) {
        int len = arg.getNumberOfPartitions();
        double logq = 0;
        double leftValue = MathUtils.nextInt(2);
        double rightValue = 1.0 - leftValue;
        logq += Math.log(2);
        if (partitioning != null) {
            int cut = MathUtils.nextInt(len - 1);
            for (int i = 0; i <= cut; i++)
                partitioning.setParameterValueQuietly(i, leftValue);
            for (int i = cut + 1; i < len; i++)
                partitioning.setParameterValueQuietly(i, rightValue);
        }
        logq += Math.log(len - 1);
        return logq;
    }


    private double AddOperation() throws ARGOperatorFailedException {

        double logq = 0;

        // Draw attachment point for new bifurcation

//	    System.err.println("Starting add operation.");

        ArrayList<NodeRef> potentialBifurcationChildren = new ArrayList<NodeRef>();
        ArrayList<NodeRef> potentialReassortmentChildren = new ArrayList<NodeRef>();

        double treeHeight = arg.getNodeHeight(arg.getRoot());

        double newBifurcationHeight = treeHeight * MathUtils.nextDouble();
        double newReassortmentHeight = treeHeight * MathUtils.nextDouble();
        if (newReassortmentHeight > newBifurcationHeight) {
            double temp = newReassortmentHeight;
            newReassortmentHeight = newBifurcationHeight;
            newBifurcationHeight = temp;
        }

//	    logq += 2.0 * Math.log(treeHeight); // Uniform before sorting

        int totalPotentialBifurcationChildren = findPotentialAttachmentPoints(
                newBifurcationHeight, potentialBifurcationChildren);
        int totalPotentialReassortmentChildren = findPotentialAttachmentPoints(
                newReassortmentHeight, potentialReassortmentChildren);

        if (totalPotentialBifurcationChildren == 0 || totalPotentialReassortmentChildren == 0) {
            throw new RuntimeException("Unable to find attachment points.");
        }

        Node recNode = (Node) potentialReassortmentChildren.get(MathUtils
                .nextInt(totalPotentialReassortmentChildren));
//	    logq += Math.log(totalPotentialReassortmentChildren);

        Node recParentL = recNode.leftParent;
        Node recParentR = recNode.rightParent;
        Node recParent = recParentL;
        if (recParentL != recParentR) {
            if (arg.getNodeHeight(recParentL) < newReassortmentHeight)
                recParent = recParentR;
            else if (arg.getNodeHeight(recParentR) > newReassortmentHeight
                    && MathUtils.nextDouble() > 0.5)
                recParent = recParentR;
        }


        Node sisNode = (Node) potentialBifurcationChildren.get(MathUtils
                .nextInt(totalPotentialBifurcationChildren));
//	    logq += Math.log(totalPotentialBifurcationChildren);

        Node sisParentL = sisNode.leftParent;
        Node sisParentR = sisNode.rightParent;
        Node sisParent = sisParentL;
        if (sisParentL != sisParentR) {
            if (arg.getNodeHeight(sisParentL) < newBifurcationHeight)
                sisParent = sisParentR;
            else if (arg.getNodeHeight(sisParentR) > newBifurcationHeight
                    && MathUtils.nextDouble() > 0.5)
                sisParent = sisParentR;
        }

        double newBifurcationRateCategory = 1.0;
        double newReassortmentRateCategory = 1.0;

        Node newBifurcation = arg.new Node();
        newBifurcation.heightParameter = new Parameter.Default(newBifurcationHeight);
        newBifurcation.rateParameter = new Parameter.Default(newBifurcationRateCategory);
        newBifurcation.setupHeightBounds();

        Node newReassortment = arg.new Node();
        newReassortment.bifurcation = false;
        newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
        newReassortment.rateParameter = new Parameter.Default(newReassortmentRateCategory);
        newReassortment.setupHeightBounds();

        arg.beginTreeEdit();
        if (sisParent.bifurcation)
            arg.singleRemoveChild(sisParent, sisNode);
        else
            arg.doubleRemoveChild(sisParent, sisNode);
        if (sisNode != recNode) {
            if (recParent.bifurcation)
                arg.singleRemoveChild(recParent, recNode);
            else
                arg.doubleRemoveChild(recParent, recNode);
        }
        if (sisParent.bifurcation)
            arg.singleAddChild(sisParent, newBifurcation);
        else
            arg.doubleAddChild(sisParent, newBifurcation);
        if (sisNode != recNode)
            arg.singleAddChild(newBifurcation, sisNode);
        arg.doubleAddChild(newReassortment, recNode);

        Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());

//        logq +=
        drawRandomPartitioning(partitioning);
//		System.err.println("point 1");

        if (sisNode != recNode) {
            arg.addChildAsRecombinant(newBifurcation, recParent,
                    newReassortment, partitioning);
        } else {
            arg.addChildAsRecombinant(newBifurcation, newBifurcation,
                    newReassortment, partitioning);
        }

        arg.expandARGWithRecombinant(newBifurcation, newReassortment,
                internalNodeParameters,
                internalAndRootNodeParameters,
                nodeRates);

//		System.err.println("point 2");

        arg.pushTreeSizeChangedEvent();

        arg.endTreeEdit();
//        try {
//            arg.checkTreeIsValid();
//        } catch (MutableTree.InvalidTreeException ite) {
//            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
//                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
//        }

//	    double d1 = findPotentialNodesToRemove(null);
//	    System.err.printf("d1 = %5.4f\n",d1);

//        logq -= Math.log(findPotentialNodesToRemove(null));
//        if (!(recNode.leftParent.isBifurcation() && recNode.rightParent.isRoot()) &&
//                !(recNode.rightParent.isBifurcation() && recNode.leftParent.isRoot()))
//            logq -= Math.log(2.0);

//	    System.err.println("End add ARG operation.");

//        int nodes = arg.getInternalNodeCount() - 1;

//
/*	    int i = findPotentialNodesToRemove(null);
	    if (i==0) {
		    System.err.println("why can't i remove this one?");
		    System.err.println("graph:"+arg.toGraphString());
		    System.exit(1);
		    return Double.NEGATIVE_INFINITY ;
	    }*/
        logq = 0;
//	    logq -= Math.log(findPotentialNodesToRemove(null));
//	    System.err.println("logq add = "+logq);
//	    System.err.println("After add:\n"+arg.toGraphString());
        return logq;
    }

//	if( !recParent1.bifurcation || recParent1.isRoot() ) { // One orientation is valid


    public void sanityCheck() {
        int len = arg.getNodeCount();
        for (int i = 0; i < len; i++) {
            Node node = (Node) arg.getNode(i);
            if (node.bifurcation) {
                boolean equalChild = (node.leftChild == node.rightChild);
                if ((equalChild && node.leftChild != null)) {
                    if (!node.leftChild.bifurcation && ((node.leftChild).leftParent == node))
                        ;
                    else {
                        System.err.println("Node " + (i + 1) + " is insane.");
                        System.err.println(arg.toGraphString());
                        System.exit(-1);
                    }
                }
            } else {
                if ((node.leftChild != node.rightChild)) {
                    System.err.println("Node " + (i + 1) + " is insane.");
                    System.err.println(arg.toGraphString());
                    System.exit(-1);
                }
            }
            if (!node.isRoot()) {
                double d;
                d = node.getHeight();
            }
        }
    }

    private int times = 0;

    private double getDelta() {
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
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

//	public int getMode() {
//		return mode;
//	}

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }


    public String getPerformanceSuggestion() {
        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(getSize(), Double.MAX_VALUE, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing size to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing size to about " + ws;
        } else return "";
    }

    public String getOperatorName() {
        return SUBTREE_SLIDE;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return SUBTREE_SLIDE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean swapRates = false;
            boolean swapTraits = false;

            double singlePartitionProbability = 0.0;
            boolean isRecombination = false;

            CoercionMode mode = CoercionMode.parseMode(xo);

//			int mode = CoercableMCMCOperator.DEFAULT;
//			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
//				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
//					mode = CoercableMCMCOperator.COERCION_ON;
//				} else {
//					mode = CoercableMCMCOperator.COERCION_OFF;
//				}
//			}

            if (xo.hasAttribute(SINGLE_PARTITION)) {
                singlePartitionProbability = xo.getDoubleAttribute(SINGLE_PARTITION);
            }

            if (xo.hasAttribute(IS_RECOMBINATION)) {
                isRecombination = xo.getBooleanAttribute(IS_RECOMBINATION);

            }

            if (xo.hasAttribute(SWAP_RATES)) {
                swapRates = xo.getBooleanAttribute(SWAP_RATES);
            }
            if (xo.hasAttribute(SWAP_TRAITS)) {
                swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
            }
            Object obj = xo.getChild(ARGModel.class);
            //TreeModel tmp = (TreeModel)xo.getChild(TreeModel.class);
            ARGModel treeModel = null;
            //if( (tmp.TREE_MODEL).compareTo(VariableSizeTreeModel.TREE_MODEL) == 0) {
            if (obj instanceof ARGModel) {
                //System.err.println("Found VSTM");
                treeModel = (ARGModel) obj;
            } else {
                System.err.println("Must specify a variable size tree model to use the AddRemoveSubtreeOperators");
                System.exit(-1);
            }

               // bug, getChild(String) returns an xmlObject;   not in XML rules
            CompoundParameter parameter1 = null; //(CompoundParameter) xo.getChild(JUST_INTERNAL);
            CompoundParameter parameter2 = null; //CompoundParameter) xo.getChild(INTERNAL_AND_ROOT);
            CompoundParameter parameter3 = null; //(CompoundParameter) xo.getChild(NODE_RATES);

            int weight = xo.getIntegerAttribute("weight");
//            int maxTips = xo.getIntegerAttribute(MAX_VALUE);
            double size = xo.getDoubleAttribute("size");
            boolean gaussian = xo.getBooleanAttribute("gaussian");
            return new ObsoleteARGNewEventOperator(treeModel, weight, size, gaussian, swapRates, swapTraits,
                    mode, parameter1, parameter2, parameter3, singlePartitionProbability, isRecombination);
        }

        public String getParserDescription() {
            return "An operator that slides a subarg.";
        }

        public Class getReturnType() {
            return ObsoleteARGAddRemoveEventOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule("weight"),
//                AttributeRule.newIntegerRule(MAX_VALUE),
                AttributeRule.newDoubleRule("size"),
                AttributeRule.newBooleanRule("gaussian"),
                AttributeRule.newBooleanRule(SWAP_RATES, true),
                AttributeRule.newBooleanRule(SWAP_TRAITS, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(ARGModel.class)//,
//			new ElementRule(Parameter.class)
        };
    };

}
