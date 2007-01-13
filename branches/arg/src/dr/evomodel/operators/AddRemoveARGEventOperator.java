/*
 * AddRemoveSubtreeOperator.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.evomodel.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ARGModel;
import dr.evomodel.tree.ARGModel.Node;
import dr.inference.model.Parameter;
import dr.inference.model.VariableSizeCompoundParameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.BitSet;

/**
 * Implements the subtree slide move.
 *
 * @author Alexei Drummond
 * @version $Id: AddRemoveARGEventOperator.java,v 1.18.2.4 2006/11/06 01:38:30 msuchard Exp $
 */
public class AddRemoveARGEventOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String SUBTREE_SLIDE = "addremoveARGEvent";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_TRAITS = "swapTraits";
    public static final String MAX_VALUE = "maxTips";

    public static final String JUST_INTERNAL = "justInternalNodes";
    public static final String INTERNAL_AND_ROOT = "internalAndRootNodes";

    private ARGModel arg = null;
    private double size = 1.0;
    private boolean gaussian = false;
    //   private boolean swapRates;
    //   private boolean swapTraits;
    private int mode = CoercableMCMCOperator.DEFAULT;
    private VariableSizeCompoundParameter internalNodeParameters;
    private VariableSizeCompoundParameter internalAndRootNodeParameters;
//	private int maxTips = 1;

    public AddRemoveARGEventOperator(ARGModel arg, int weight, double size, boolean gaussian,
                                     boolean swapRates, boolean swapTraits, int mode, int maxTips,
                                     VariableSizeCompoundParameter param1,
                                     VariableSizeCompoundParameter param2) {
        this.arg = arg;
        setWeight(weight);
//		this.maxTips = maxTips;
        this.size = size;
        this.gaussian = gaussian;
//        this.swapRates = swapRates;
//        this.swapTraits = swapTraits;
        this.internalNodeParameters = param1;
        this.internalAndRootNodeParameters = param2;

        this.mode = mode;
    }


    /**
     * Do a add/remove reassortment node operation
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() throws OperatorFailedException {
//		System.err.println("Starting AddRemove Operation");
        if (MathUtils.nextDouble() < 0.5)
            return AddOperation();
        else
            return RemoveOperation();
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

    private double RemoveOperation() throws OperatorFailedException {
        double logq = 0;

        // 1. Draw reassortment node uniform randomly

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();

        int totalPotentials = findCurrentReassortmentNodes(potentialNodes);
        if (totalPotentials == 0)
            throw new OperatorFailedException("No reassortment nodes to remove.");
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
                try {
                    arg.endTreeEdit();
                } catch (MutableTree.InvalidTreeException ite) {
                    throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                            + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
                }
                throw new OperatorFailedException("Not reversible deletion.");
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
            // TODO should be able to remove at some point
            if ((recChild.getHeight() > recParent2.getHeight()) ||
                    (otherChild.getHeight() > recGrandParent.getHeight())) {
                try {
                    arg.endTreeEdit();
                } catch (MutableTree.InvalidTreeException ite) {
                    //throw new OperatorFailedException("Not reversible deletion.");
                }
                System.err.println("How did I get here?");
                System.exit(-1);
            }
            recParent = recParent1;
        }
        if (doneSomething) {
//			System.err.println("Trying to remove "+recNode.number+" and "+recParent.number);
//		
            arg.contractARGWithRecombinant(recNode, recParent, internalNodeParameters, internalAndRootNodeParameters);
        }
//		System.err.println("End ARG\n"+arg.toGraphString());
        arg.pushTreeSizeChangedEvent();
        try {
            arg.endTreeEdit();
        } catch (MutableTree.InvalidTreeException ite) {
            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
        }
        logq -= Math.log(reverseBifurcationSpan) - Math.log(reverseReassortmentSpan);
        logq -= Math.log(arg.getNodeCount() + arg.getReassortmentNodeCount() - 1);
        logq -= Math.log(this.findPotentialAttachmentSisters(recChild,
                arg.getNodeHeight(recChild), null));
        return -logq; // 1 / total potentials * 1 / 2 (if valid) * length1 * length2 * attachmentSisters
    }

    private void checkAllHeights() {
        int len = arg.getInternalNodeCount();
        System.err.println("# internal nodes = " + len);
        int n = internalNodeParameters.getNumParameters();
        System.err.println("VSCP (" + n + ")");
        for (int i = 0; i < n; i++) {
            System.err.println(internalNodeParameters.getParameterValue(i));
        }
        n = arg.getInternalNodeCount();
        System.err.println("Checking all internal nodes (" + n + ") via tree:");
        for (int i = 0; i < n; i++) {
            NodeRef node = arg.getInternalNode(i);
            System.err.print(Tree.Utils.uniqueNewick(arg, node) + " ");
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

    private int drawOnePartition(BitSet partitionSet) {
        int draw = MathUtils.nextInt(partitionSet.cardinality());
//		System.err.println("draw = "+draw);
        int result = partitionSet.nextSetBit(0);
        for (int i = 0; i < draw; i++)
            result = partitionSet.nextSetBit(result + 1);
        return result;
    }

    private BitSet bsTransportor = null;
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
//				// [TODO -- allow rerooting
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


    private double AddOperation() throws OperatorFailedException {
        double logq = 0;

        // Draw potential places to add a reassortment node

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
        int totalPotentials = findPotentialReassortmentNodes(potentialNodes);
        if (totalPotentials == 0)
            throw new OperatorFailedException("No more nodes to recombine.");
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
        if (totalAttachments == 0)
            throw new OperatorFailedException("no more attachment points for this recomb");
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

        BitSet bitLeft = new BitSet();
        BitSet bitRight = new BitSet();
        if (MathUtils.nextDouble() < 0.5) {
            bitLeft.set(0);
            bitRight.set(1);
        } else {
            bitLeft.set(1);
            bitRight.set(0);
        }
        logq += Math.log(2.0);

        if (sisNode != recNode) {
            arg.addChildAsRecombinant(newBifurcation, recParent,
                    newReassortment, bitLeft, bitRight);
        } else {
            arg.addChildAsRecombinant(newBifurcation, newBifurcation,
                    newReassortment, bitLeft, bitRight);
        }
        arg.expandARGWithRecombinant(newBifurcation, newReassortment,
                internalNodeParameters,
                internalAndRootNodeParameters);
        // arg.addNewHeightParameter(newReassortment.heightParameter,internalNodeParameters);
        // arg.expandNodesWithRecombinant(arg.getRoot(), recNR);
        // arg.expandNodesWithSubtree(arg.getRoot(), newSubtree);
        // arg.reconstructTrees();
//		System.err.println("recNode   = " + recNode.number);
//		System.err.println("recParent = " + recParent.number);
//		System.err.println("sisNode   = " + sisNode.number);
//		System.err.println("sisParent = " + sisParent.number);
//		System.err.println("End ARG\n" + arg.toGraphString());
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
        try {
            arg.endTreeEdit();
        } catch (MutableTree.InvalidTreeException ite) {
            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                    + "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
        }
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

        return -logq;
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

    public int getMode() {
        return mode;
    }

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

            int mode = CoercableMCMCOperator.DEFAULT;
            if (xo.hasAttribute(AUTO_OPTIMIZE)) {
                if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
                    mode = CoercableMCMCOperator.COERCION_ON;
                } else {
                    mode = CoercableMCMCOperator.COERCION_OFF;
                }
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
//        	VariableSizeCompoundParameter parameter = (VariableSizeCompoundParameter)xo.getChild(Parameter.class);

            VariableSizeCompoundParameter parameter1 = (VariableSizeCompoundParameter) xo.getSocketChild(JUST_INTERNAL);
            VariableSizeCompoundParameter parameter2 = (VariableSizeCompoundParameter) xo.getSocketChild(INTERNAL_AND_ROOT);

            //VariableSizeTreeModel treeModel = (VariableSizeTreeModel)xo.getChild(TreeModel.class);
            int weight = xo.getIntegerAttribute("weight");
            int maxTips = xo.getIntegerAttribute(MAX_VALUE);
            double size = xo.getDoubleAttribute("size");
            boolean gaussian = xo.getBooleanAttribute("gaussian");
            return new AddRemoveARGEventOperator(treeModel, weight, size, gaussian, swapRates, swapTraits,
                    mode, maxTips, parameter1, parameter2);
        }

        public String getParserDescription() {
            return "An operator that slides a subarg.";
        }

        public Class getReturnType() {
            return AddRemoveARGEventOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule("weight"),
                AttributeRule.newIntegerRule(MAX_VALUE),
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
