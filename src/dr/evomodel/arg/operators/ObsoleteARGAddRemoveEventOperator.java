/*
 * ObsoleteARGAddRemoveEventOperator.java
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
import dr.evolution.tree.Tree;
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
public class ObsoleteARGAddRemoveEventOperator extends AbstractCoercableOperator {
//		SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String SUBTREE_SLIDE = "addremoveARGEvent";
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

    public ObsoleteARGAddRemoveEventOperator(ARGModel arg, int weight, double size, boolean gaussian,
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
        } catch (ARGOperatorFailedException ofe) {
            if (ofe.getMessage().compareTo("No reassortment nodes to remove.") != 0) {
                System.err.println("Catch: " + ofe.getMessage());
//            System.exit(-1);
            }
        }
        if (arg.isBifurcationDoublyLinked(arg.getRoot()))
            throw new RuntimeException("trouble with double-rooted root");
        return logq;
    }

/*
    private ArrayList<NodeRef> getPotentialSubtreesToMerge() {
        ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
        int n = arg.getNodeCount();
        for(int i=0; i<n; i++) {
            Node node = (Node)arg.getNode(i);
            if( node.dupSister != null ) {
                potentials.add(node);
            }
        }
        return potentials;
    }

    private int countPotentialSubtreesToMerge() {
        int count = 0;
        //ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
        int n = arg.getNodeCount();
        for(int i=0; i<n; i++) {
            Node node = (Node)arg.getNode(i);
            if( node.dupSister != null ) {
                //potentials.add(node);
                count++;
            }
        }
        //return potentials;
        return count;
    }
    */

    /*
                * Generates a list of all branches in the ARG, there each branch leads up the ARG from the saved nodes
                *
                * Should return a deterministic function of # bifurcations and # reassortments
                *
                */

    private int findPotentialReassortmentNodes(ArrayList<NodeRef> list) {
        int count = 0;
        int n = arg.getNodeCount();
        for (int i = 0; i < n; i++) {
            Node node = (Node) arg.getNode(i);
            if (node.leftParent != null) { // rules out only the root; may subclass for different types of problems
                list.add(node);
                count++;
                if (node.isReassortment()) {
                    list.add(node);
                    count++; // reassortment nodes have two parent branches
                }
            }
        }
//		int check = arg.getNodeCount() + arg.getReassortmentNodeCount() - 1;
//		if( check != count ) {
//			System.err.println("What the fuck?");
//			System.exit(-1);
//		}
        return count;
    }

    private int findCurrentReassortmentNodes(ArrayList<NodeRef> list) {
        int count = 0;
        int n = arg.getNodeCount();
        Node root = (Node) arg.getRoot();
        for (int i = 0; i < n; i++) {
            Node node = (Node) arg.getNode(i);
            if (node.isReassortment() && (node.leftParent != root && node.rightParent != root)) {
                if (list != null)
                    list.add(node);
                count++;
            }
        }
        return count;
    }

    /*private double RemoveOperation() throws OperatorFailedException {
             arg.beginTreeEdit();
             if (arg.getReassortmentNodeCount() > 0)
                 arg.removeNullCounter();
              try {
                         arg.endTreeEdit();
                     } catch (MutableTree.InvalidTreeException ite) {
                         throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                                 + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
                     }
             arg.pushTreeSizeChangedEvent();
             return 0;
         }*/

    /*   private double AddOperation() throws OperatorFailedException {
            arg.beginTreeEdit();
            arg.addNullCounter();
            try {
                        arg.endTreeEdit();
                    } catch (MutableTree.InvalidTreeException ite) {
                        throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                                + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
                    }
            arg.pushTreeSizeChangedEvent();
            return 0;
        }
    */

    private double RemoveOperation() throws ARGOperatorFailedException {
        double logq = 0;

//	    System.err.println("Starting remove ARG operation.");

        // 1. Draw reassortment node uniform randomly

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();

        int totalPotentials = findCurrentReassortmentNodes(potentialNodes);
        if (totalPotentials == 0)
            throw new ARGOperatorFailedException("No reassortment nodes to remove.");
        Node recNode = (Node) potentialNodes.get(MathUtils.nextInt(totalPotentials));
        logq += Math.log(totalPotentials);

        double reverseReassortmentSpan = 0;
        double reverseBifurcationSpan = 0;

        arg.beginTreeEdit();
        boolean doneSomething = false;
        Node recParent = recNode.leftParent;
        Node recChild = recNode.leftChild;
        if (recNode.leftParent == recNode.rightParent) { // Doubly linked.
            Node recGrandParent = recParent.leftParent;

            reverseReassortmentSpan = arg.getNodeHeight(recParent) - arg.getNodeHeight(recChild);
            reverseBifurcationSpan = arg.getNodeHeight(recGrandParent) - arg.getNodeHeight(recChild);
            if (arg.isRoot(recParent)) { // This case should never happen as double links
                arg.setRoot(recChild);    // to root can not be added or removed.
            } else {
                //Node recGrandParent = recParent1.leftParent; // And currently recParent must be a bifurcatio
//				System.err.println("recGrand   ="+recGrandParent.number);
                arg.doubleRemoveChild(recGrandParent, recParent);
                arg.doubleRemoveChild(recNode, recChild);
                if (recGrandParent.bifurcation)
                    arg.singleAddChild(recGrandParent, recChild);
                else
                    arg.doubleAddChild(recGrandParent, recChild);
            }
            doneSomething = true;
            // There are not left/right choices to be made for doubly linked removals.
            // End doubly linked.
        } else { // Two different parents.
            Node recParent1 = recNode.leftParent;
            Node recParent2 = recNode.rightParent;
            if ((!recParent1.bifurcation && !recParent2.bifurcation) ||
                    (!recParent1.bifurcation && recParent2.isRoot()) ||
                    (!recParent2.bifurcation && recParent1.isRoot())) {

                arg.endTreeEdit();

//                try {
//                    arg.checkTreeIsValid();
//                } catch (MutableTree.InvalidTreeException ite) {
//                    throw new RuntimeException(ite.toString() + "\n" + arg.toString()
//                            + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
//                }
                throw new ARGOperatorFailedException("Not reversible deletion.");
            }
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
//			System.err.println("recNode    ="+recNode.number);
//			System.err.println("recParent1 ="+recParent1.number);
//			System.err.println("recParent2 ="+recParent2.number);

            Node recGrandParent = recParent1.leftParent; // And currently recParent must be a bifurcatio
//			System.err.println("recGrand   ="+recGrandParent.number);
            //Node recChild = recNode.leftChild;
            Node otherChild = recParent1.leftChild;
            if (otherChild == recNode)
                otherChild = recParent1.rightChild;
//			System.err.println("recChild   ="+recChild.number);
//			System.err.println("otherChild ="+otherChild.number);
            reverseReassortmentSpan = Math.min(arg.getNodeHeight(recParent1), arg.getNodeHeight(recParent2))
                    - arg.getNodeHeight(recChild);
            reverseBifurcationSpan = arg.getNodeHeight(recGrandParent)
                    - Math.max(arg.getNodeHeight(recChild), arg.getNodeHeight(otherChild));
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
            // Check for node height troubles

            if ((recChild.getHeight() > recParent2.getHeight()) ||
                    (otherChild.getHeight() > recGrandParent.getHeight())) {
                arg.endTreeEdit();
                throw new RuntimeException("How did I get here?");
            }
            recParent = recParent1;
        }
        if (doneSomething) {
//			System.err.println("Trying to remove "+recNode.number+" and "+recParent.number);
//
            arg.contractARGWithRecombinant(recParent, recNode,
                    internalNodeParameters, internalAndRootNodeParameters, nodeRates);
        }
//		System.err.println("End ARG\n"+arg.toGraphString());


        arg.pushTreeSizeChangedEvent();
        arg.endTreeEdit();
//        try {
//            arg.checkTreeIsValid();
//        } catch (MutableTree.InvalidTreeException ite) {
//            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
//                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
//        }

//	    System.err.println("Checking remove validity.");

        // todo -- check all ARGTree.Roots
//	    if (!arg.validRoot())
//	        throw new OperatorFailedException("Roots are invalid");

//         logq -= Math.log(reverseBifurcationSpan) + Math.log(reverseReassortmentSpan);   // TODO removed? because of prior ratio
        logq -= Math.log(arg.getNodeCount() + arg.getReassortmentNodeCount() - 1); // findPotentialRessortmentNodes()
        logq -= Math.log(this.findPotentialAttachmentSisters(recChild,
                arg.getNodeHeight(recChild), null));
//        System.err.println(drawRandomPartitioning(null));
        // logq -= drawRandomPartitioning(null);

//	    System.err.println("End remove ARG operation.");


        int nodes = arg.getInternalNodeCount() - 1;
        logq += lnGamma(nodes) - lnGamma(nodes + 2); // TODO move into prior

        logq += 3 * Math.log(2);

        logq = 0;

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

    /*private void checkAllPartitionLabels() {
                   int len = arg.getExternalNodeCount();
                   System.err.println("Tip Partitions:");
                   for(int i=0; i<len; i++) {
                       Node node = (Node)arg.getExternalNode(i);
                       BitSet partSet = node.partitionSet;
                       System.err.print(Tree.Utils.uniqueNewick(arg, node)+":");
                       for(int j=partSet.nextSetBit(0); j>=0; j=partSet.nextSetBit(j+1)) {
                           System.err.print(" "+j);
                       }
                       System.err.println();
                   }
               }*/

/*	private ArrayList<NodeRef> getPotentialTipsToMerge() {
        ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
        int len = arg.getExternalNodeCount();
        for(int i=0; i<len; i++){
            NodeRef nr = arg.getExternalNode(i);
            BitSet pS = ((Node)nr).partitionSet;
            if( pS.cardinality() == 1)
                potentials.add(nr);
        }
        return potentials;
    }

    private ArrayList<Node> getSameTaxonTips(Node asMe) {
        ArrayList<Node> potentials = new ArrayList<Node>();
        int len = arg.getExternalNodeCount();
        Taxon same = asMe.taxon;
        for(int i=0; i<len; i++){
            Node n = (Node)arg.getExternalNode(i);
            if( (n != asMe) && (n.taxon == same) )
                potentials.add(n);
        }
        return potentials;
    }

    private ArrayList<NodeRef> getPotentialTipsToSplit() {
        ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
        int len = arg.getExternalNodeCount();
        for(int i=0; i<len; i++) {
            NodeRef nr = arg.getExternalNode(i);
            BitSet pS = ((Node)nr).partitionSet;
            if( pS.cardinality() > 1 )  // can be split
                potentials.add(nr);
        }
        return potentials;
    }*/

    /*;private int drawOnePartition(BitSet partitionSet) {
             int draw = MathUtils.nextInt(partitionSet.cardinality());
     //		System.err.println("draw = "+draw);
             int result = partitionSet.nextSetBit(0);
             for (int i = 0; i < draw; i++)
                 result = partitionSet.nextSetBit(result + 1);
             return result;
         }*/

//    private BitSet bsTransportor = null;
    /*
               private ArrayList<NodeRef> getPotentialSubtreesToSplit() {
                   ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
                   // Start with all possible nodes except the root
                   BitSet bsTmp = null;
                   int len = arg.getNodeCount();
                   for(int i=0; i<len; i++) {
                       NodeRef nr = arg.getNode(i);
                       if( ! arg.isRoot(nr) ) {
                           // All tips of nr must have at least two partitions
                           ArrayList<NodeRef> allTips = arg.getDescendantTipNodes(nr);
                           int n = allTips.size();
                           boolean valid = true;
                           for(int j=0; valid && j<n; j++) {
                               BitSet pS = ((Node)allTips.get(j)).partitionSet;
                               if( j == 0 ) {
                                   bsTmp = (BitSet)pS.clone();
                               }
                               else
                                   bsTmp.and(pS);
                               if( pS.cardinality() == 1 )  // can be split
                                   valid = false;
                           }
                           if( valid && (bsTmp.cardinality()>1) )
                               potentials.add(nr);
                       }
                   }
                   bsTransportor = bsTmp;
                   return potentials;
               }
               */
/*	private int countPotentialSubtreesToSplit() {
//		ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
		int count = 0;
		// Start with all possible nodes except the root
		BitSet bsTmp = null;
		int len = arg.getNodeCount();
		for(int i=0; i<len; i++) {
			NodeRef nr = arg.getNode(i);
			if( ! arg.isRoot(nr) ) {
				// All tips of nr must have at least two partitions
				ArrayList<NodeRef> allTips = arg.getDescendantTipNodes(nr);
				int n = allTips.size();
				boolean valid = true;
				for(int j=0; valid && j<n; j++) {
					BitSet pS = ((Node)allTips.get(j)).partitionSet;
					if( j == 0 ) {
						bsTmp = (BitSet)pS.clone();
					}
					else
						bsTmp.and(pS);
					if( pS.cardinality() == 1 )  // can be split
						valid = false;
				}
				if( valid && (bsTmp.cardinality()>1) )
					//potentials.add(nr);
					count++;
			}
		}
		bsTransportor = bsTmp;
		//return potentials;
		return count;
	}*/

//	private ArrayList<NodeRef> getPotentialReattachments(double min) {
//		ArrayList<NodeRef> potentials = new ArrayList<NodeRef>();
//		int len = arg.getNodeCount();
//		for(int i=0; i<len; i++) {
//			NodeRef nr = arg.getNode(i);
//			if( !arg.isRoot(nr) && (arg.getMinParentNodeHeight(nr) > min) ) {
//				// Do not reroot and reattach only above current height
//
//				potentials.add(nr);
//			}
//		}
//		return potentials;
//	}

    private int findPotentialAttachmentSisters(NodeRef rec, double min, ArrayList<NodeRef> list) {
        int count = 0;
        int len = arg.getNodeCount();
        for (int i = 0; i < len; i++) {
            Node nr = (Node) arg.getNode(i);
            if (!nr.isRoot()) {
                if (arg.getNodeHeight(nr.leftParent) > min) {
                    // can add node somewhere between min and nr.parent
                    if (list != null)
                        list.add(nr);
                    count++;
                }
                if (nr.isReassortment() && arg.getNodeHeight(nr.rightParent) > min) {
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
            boolean first = MathUtils.nextBoolean();
            if (partitioning != null) {
                if (first)
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

//	    System.err.println("Starting add ARG operation. ... "+arg.getReassortmentNodeCount());
        // Draw potential places to add a reassortment node

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
        int totalPotentials = findPotentialReassortmentNodes(potentialNodes);
        if (totalPotentials == 0) {
            System.err.println("Should never get here AA");
            System.exit(-1);
            throw new ARGOperatorFailedException("No more nodes to recombine.");
        }
        Node recNode = (Node) potentialNodes.get(MathUtils
                .nextInt(totalPotentials));
        Node recParentL = recNode.leftParent;
        Node recParentR = recNode.rightParent;
        Node recParent = recParentL;
        if ((recParentL != recParentR) && MathUtils.nextDouble() > 0.5)
            recParent = recParentR;

        logq += Math.log(totalPotentials);

        double minHeight = arg.getNodeHeight(recNode); // Attachment must occur
        // previous to this time

        ArrayList<NodeRef> attachments = new ArrayList<NodeRef>();
        int totalAttachments = findPotentialAttachmentSisters(recNode,
                minHeight, attachments);
        if (totalAttachments == 0) {
//            System.err.println("Should never get here AB");
//            System.exit(-1);

            throw new ARGOperatorFailedException("no more attachment points for this recomb");
        }
        Node sisNode = (Node) attachments.get(MathUtils
                .nextInt(totalAttachments));
        Node sisParentL = sisNode.leftParent;
        Node sisParentR = sisNode.rightParent;
        Node sisParent = sisParentL;
        if (sisParentL != sisParentR) {
            if (arg.getNodeHeight(sisParentL) <= minHeight)
                sisParent = sisParentR;
            else if (arg.getNodeHeight(sisParentR) > minHeight
                    && MathUtils.nextDouble() > 0.5)
                sisParent = sisParentR;
        }

        logq += Math.log(totalAttachments);

        Node newBifurcation = arg.new Node();
        Node newReassortment = arg.new Node();
        newReassortment.bifurcation = false;
        double sisHeight = sisNode.getHeight();
        if (sisHeight < minHeight) {
            sisHeight = minHeight;
        }
        double spHeight = arg.getNodeHeight(sisParent);
        double totalLength = spHeight - sisHeight;

        double newLength = sisHeight + MathUtils.nextDouble() * totalLength;

        logq -= Math.log(totalLength); // prior ratio

        newBifurcation.heightParameter = new Parameter.Default(newLength);
        newBifurcation.setupHeightBounds();

        logq += Math.log(totalLength); // Uniform[spHeight, sisHeight]

        double topHeight = newLength;
        double recParentHeight = arg.getNodeHeight(recParent);
        if (topHeight > recParentHeight)
            topHeight = recParentHeight;
        double recHeight = arg.getNodeHeight(recNode);
        totalLength = topHeight - recHeight;
        newLength = recHeight + MathUtils.nextDouble() * totalLength;

        logq -= Math.log(totalLength); // prior ratio

        newReassortment.heightParameter = new Parameter.Default(newLength);
        newReassortment.setupHeightBounds();
        logq += Math.log(totalLength); // Uniform[bifurcationHeight,recNodeHeight]

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

/*		BitSet bitLeft = new BitSet();
		BitSet bitRight = new BitSet();
		if (MathUtils.nextDouble() < 0.5) {
			bitLeft.set(0);
			bitRight.set(1);
		} else {
			bitLeft.set(1);
			bitRight.set(0);
		}
		logq += Math.log(2.0);*/
        Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());

//        logq += drawRandomPartitioning(partitioning);          // TODO add back in

        /*    System.err.print("p: ");
                  double[] v = partitioning.getParameterValues();
                  for (double d : v)
                      System.err.print(d+" ");
                  System.err.println("");*/


        if (sisNode != recNode) {
            arg.addChildAsRecombinant(newBifurcation, recParent,
                    newReassortment, partitioning);
        } else {
            arg.addChildAsRecombinant(newBifurcation, newBifurcation,
                    newReassortment, partitioning);
        }

        /* if (sisNode != recNode) {
                        arg.addChildAsRecombinant(newBifurcation, recParent,
                                newReassortment, bitLeft, bitRight);
                    } else {
                        arg.addChildAsRecombinant(newBifurcation, newBifurcation,
                                newReassortment, bitLeft, bitRight);
                    }*/
        arg.expandARGWithRecombinant(newBifurcation, newReassortment,
                internalNodeParameters,
                internalAndRootNodeParameters,
                nodeRates);
        // arg.addNewHeightParameter(newReassortment.heightParameter,internalNodeParameters);
        // arg.expandNodesWithRecombinant(arg.getRoot(), recNR);
        // arg.expandNodesWithSubtree(arg.getRoot(), newSubtree);
        // arg.reconstructTrees();
//
//		System.err.println("End ARG\n" + arg.toGraphString()); System.err.println("recNode   = " + recNode.number);
//		System.err.println("recParent = " + recParent.number);
//		System.err.println("sisNode   = " + sisNode.number);
//		System.err.println("sisParent = " + sisParent.number);
        //System.exit(-1);
//		ARGTree tree = new ARGTree(arg);
        // System.err.println("Reassortment nodes =
        // "+arg.getReassortmentNodeCount());
        // //System.exit(-1);
        // //System.err.println("ARGTree = "+Tree.Utils.uniqueNewick(tree,
        // tree.getRoot()));
//		tree = new ARGTree(arg, 0);
//		System.err.println("Part 0  = " + tree.toString());
//		tree = new ARGTree(arg, 1);
//		System.err.println("Part 1  = " + tree.toString());
        // //if( recNode == sisNode)
        //System.exit(-1);
//    	System.err.println("Sanity check in Add Operation");
//		sanityCheck();
//		System.err.println("End Add Operator sanity check");

        arg.pushTreeSizeChangedEvent();

        arg.endTreeEdit();

//        try {
//            arg.checkTreeIsValid();
//        } catch (MutableTree.InvalidTreeException ite) {
//            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
//                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
//        }

//	    System.err.println("Checking add validity.");
        // todo -- check all ARGTree.Roots
//	    if (!arg.validRoot())
//	        throw new OperatorFailedException("Roots are invalid");

        // System.err.println("Made it thru once!");

        // System.exit(-1);
        // times++;
        // System.err.println("Adds = "+times);
        // if( times >= 4 )
        // throw new OperatorFailedException("Do many tries!");
        // System.err.println("Add a recomb node.");
//		int cnt = arg.getReassortmentNodeCount();
//		if (cnt > 20)
//			throw new OperatorFailedException("No more than X reassortments");

        logq -= Math.log(findCurrentReassortmentNodes(null));
        if (!(recNode.leftParent.isBifurcation() && recNode.rightParent.isRoot()) &&
                !(recNode.rightParent.isBifurcation() && recNode.leftParent.isRoot()))
            logq -= Math.log(2.0);

//	    System.err.println("End add ARG operation.");

        int nodes = arg.getInternalNodeCount() - 1;
        logq += lnGamma(nodes) - lnGamma(nodes - 2); // TODO move into prior

        logq = 0;
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

    private int intersectingEdges(Tree tree, NodeRef node, double height, ArrayList directChildren) {

        NodeRef parent = arg.getParent(node);

        if (arg.getNodeHeight(parent) < height) return 0;

        if (arg.getNodeHeight(node) < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < arg.getChildCount(node); i++) {
            count += intersectingEdges(tree, arg.getChild(node, i), height, directChildren);
        }
        return count;
    }

    /**
     * @return the other child of the given parent.
     */
    private NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

        if (arg.getChild(parent, 0) == child) {
            return arg.getChild(parent, 1);
        } else {
            return arg.getChild(parent, 0);
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

//			int mode = CoercableMCMCOperator.DEFAULT;
//			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
//				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
//					mode = CoercableMCMCOperator.COERCION_ON;
//				} else {
//					mode = CoercableMCMCOperator.COERCION_OFF;
//				}
//			}

            CoercionMode mode = CoercionMode.parseMode(xo);

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
            return new ObsoleteARGAddRemoveEventOperator(treeModel, weight, size, gaussian, swapRates, swapTraits,
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
