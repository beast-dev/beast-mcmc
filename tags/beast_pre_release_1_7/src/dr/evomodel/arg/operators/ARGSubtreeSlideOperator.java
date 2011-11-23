/*
 * RGSubtreeSlideOperator.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.evomodel.arg.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Implements the subtree slide move.
 *
 * @author Alexei Drummond
 * @version $Id: ARGSubtreeSlideOperator.java,v 1.1.2.2 2006/11/06 01:38:30 msuchard Exp $
 */
public class ARGSubtreeSlideOperator extends AbstractCoercableOperator {
//		SimpleMCMCOperator implements CoercableMCMCOperator {

	public static final String SUBTREE_SLIDE = "argSubtreeSlide";
	//	public static final String
	public static final String SWAP_RATES = "swapRates";
	public static final String SWAP_TRAITS = "swapTraits";
	public static final String DIRICHLET_BRANCHES = "branchesAreScaledDirichlet";

	private ARGModel tree = null;
	private double size = 1.0;
	private boolean gaussian = false;
	private boolean swapRates;
	private boolean swapTraits;
	private boolean scaledDirichletBranches;
//	private int mode = CoercableMCMCOperator.DEFAULT;
//	CoercionMode model;

	public ARGSubtreeSlideOperator(ARGModel tree, int weight, double size, boolean gaussian, boolean swapRates,
	                               boolean swapTraits, boolean scaledDirichletBranches, CoercionMode mode) {
		super(mode);
		this.tree = tree;
		setWeight(weight);

		this.size = size;
		this.gaussian = gaussian;
		this.swapRates = swapRates;
		this.swapTraits = swapTraits;
		this.scaledDirichletBranches = scaledDirichletBranches;

//		this.mode = mode;
	}

	public void sanityCheck() {
		int len = tree.getNodeCount();
		for (int i = 0; i < len; i++) {
			Node node = (Node) tree.getNode(i);
			if (node.bifurcation) {
				boolean equalChild = (node.leftChild == node.rightChild);
				if ((equalChild && node.leftChild != null)) {
					if (!node.leftChild.bifurcation && ((node.leftChild).leftParent == node))
						;
					else {
						System.err.println("Node " + (i + 1) + " is insane.");
						System.err.println(tree.toGraphString());
						System.exit(-1);
					}
				}
			} else {
				if ((node.leftChild != node.rightChild)) {
					System.err.println("Node " + (i + 1) + " is insane.");
					System.err.println(tree.toGraphString());
					System.exit(-1);
				}
			}
		}
	}

	/**
	 * Do a probablistic subtree slide move.
	 *
	 * @return the log-transformed hastings ratio
	 */
	public double doOperation() throws OperatorFailedException {

//		System.err.println("Starting Subtree Slide Operation.");
		double logq = 0;

		double oldTreeHeight = tree.getNodeHeight(tree.getRoot());

		NodeRef i, newParent, newChild;

		// 1. choose a random node avoiding root
		ArrayList<NodeRef> potentialSubtrees = new ArrayList<NodeRef>();
		int numPotentialSubtrees = this.getSlideableSubtrees(tree, potentialSubtrees);
		//      	System.err.println("Slide:\n"+tree.toGraphString());
		i = potentialSubtrees.get(MathUtils.nextInt(numPotentialSubtrees));

		//       logq = - Math.log(numPotentialSubtrees);

		NodeRef iP = tree.getParent(i);

		// TODO Start rewriting here.
		NodeRef CiP = getOtherChild(tree, iP, i);
		NodeRef PiP;
		if (tree.isBifurcation(iP))
			PiP = tree.getParent(iP);
		else {
			PiP = tree.getParent(iP, MathUtils.nextInt(2));
			logq -= Math.log(2); // TODO check if really necessary.
		}

		// 2. choose a delta to move
		double delta = getDelta();
		double oldHeight = tree.getNodeHeight(iP);
		double newHeight = oldHeight + delta;

		//newHeight = tree.getNodeHeight(tree.getRoot()) + delta;

		// 3. if the move is up
		if (delta > 0) {

			// 3.1 if the topology will change
			if (PiP != null && tree.getNodeHeight(PiP) < newHeight) {

				// find new parent
				newParent = PiP;
				newChild = iP;
				while (tree.getNodeHeight(newParent) < newHeight) {
					newChild = newParent;
					if (tree.isBifurcation(newParent))
						newParent = tree.getParent(newParent);
					else {
						newParent = tree.getParent(newParent, MathUtils.nextInt(2));
						logq -= Math.log(2); // TODO check if correct.
					}
					if (newParent == null) break;
				}

				//System.err.println("No problem climbing");

//                logq += Math.log(numPotentialSubtrees);

				tree.beginTreeEdit();

				// 3.1.1 if creating a new root
				if (tree.isRoot(newChild)) {

					if (true) {
                        tree.endTreeEdit();
						try {
                            tree.checkTreeIsValid();
						} catch (MutableTree.InvalidTreeException e) {
							e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
						}
						throw new OperatorFailedException("Temporarily disable re-rooting");
					}

					//Parameter rootParameter = ((Node)newChild).heightParameter;
					//Parameter otherParameter = ((Node)iP).heightParameter;
					//tree.swapHeightParameters(newChild,iP);
					tree.doubleRemoveChild(iP, CiP); // iP is always bifurcation
					tree.doubleRemoveChild(PiP, iP); // PiP can be reassortment
					tree.doubleAddChild(iP, newChild); // iP and newChild are always bifurcations, i still connected to iP
					if (tree.isBifurcation(PiP))
						tree.addChild(PiP, CiP);
					else
						tree.addChild(PiP, CiP);
					tree.setRoot(iP);
					//System.err.println("Creating new root!");
				}
				// 3.1.2 no new root
				else {
					boolean doubly = tree.isBifurcationDoublyLinked(newParent);

					tree.doubleRemoveChild(iP, CiP); // iP bifurcation
					tree.doubleRemoveChild(PiP, iP); // PiP can be reassortment
					tree.doubleRemoveChild(newParent, newChild);
					//tree.doubleAddChild(iP, newChild);

					tree.doubleAddChild(PiP, CiP);
					if (!doubly) {
						tree.doubleAddChild(iP, newChild);
						tree.doubleAddChild(newParent, iP);
					} else {
						tree.singleAddChild(newParent, iP);
						tree.singleAddChildWithOneParent(iP, newChild);
						tree.singleAddChild(newParent, newChild);
					}
					//                   System.err.println("No new root!");
				}

				//               System.err.println("i  ="+((Node)i).number);
				//               System.err.println("iP  ="+((Node)iP).number);
				//               if( newParent != null )
				//               System.err.println("newParent ="+((Node)newParent).number+" "+tree.getNodeHeight(newParent));
				//               else
				//               	System.err.println("newParent is above root");
				//               System.err.println("newChild  ="+((Node)newChild).number+" "+tree.getNodeHeight(newChild));


				tree.setNodeHeight(iP, newHeight);
				//               System.err.println("iP height ="+tree.getNodeHeight(iP));

				//               System.err.println("Intermediate slide up:\n"+tree.toGraphString());

                tree.endTreeEdit();
				try {
                    tree.checkTreeIsValid();
				} catch (MutableTree.InvalidTreeException ite) {
					throw new RuntimeException(ite.toString());
				}

				// 3.1.3 count the hypothetical sources of this destination.
				int possibleSources = intersectingEdges(tree, newChild, iP, oldHeight, null);
				//               System.err.println("possible sources = " + possibleSources);

				logq -= Math.log(possibleSources);
			} else {
				// 3.2
				// just change the node height
				tree.setNodeHeight(iP, newHeight);
				logq = 0.0; // TODO check is losing -Log(2) is correct.
			}

//			System.err.println("Sanity check up-slide");
//			sanityCheck();
		}
		// 4 if we are sliding the subtree down.
		else {
			//          logq = 0;

			// 4.0 is it a valid move?
			if (tree.getNodeHeight(i) > newHeight) {
				return Double.NEGATIVE_INFINITY;
			}

			// 4.1 will the move change the topology
			if (tree.getNodeHeight(CiP) > newHeight) {
				//           	System.err.println("Starting down-slide:\n"+tree.toGraphString());
				ArrayList<NodeRef[]> newChildren = new ArrayList<NodeRef[]>();
				int possibleDestinations = intersectingEdges(tree, CiP, iP, newHeight, newChildren);

				// if no valid destinations then return a failure
				if (newChildren.size() == 0) //{ return Double.NEGATIVE_INFINITY; }
					throw new OperatorFailedException("no valid destinations");

				// pick a random parent/child destination edge uniformly from options
				int childIndex = MathUtils.nextInt(newChildren.size());
				NodeRef[] draw = (NodeRef[]) newChildren.get(childIndex);
				newChild = draw[1];
				newParent = draw[0];
				//	((NodeRef[])newChildren.get(childIndex))[1];
				//int choice = MathUtils.nextInt(2);
				NodeRef oops = null;
//                if( tree.isBifurcation(newChild) )
//                	newParent = tree.getParent(newChild);
//                else {
//                	//newParent = tree.getParent(newChild,MathUtils.nextInt(2));
//                	newParent = tree.getParent(newChild,choice);
//                	if( tree.getNodeHeight(newParent) < newHeight )
//                		newParent = tree.getParent(newChild,1-choice);
//                	oops = tree.getParent(newChild,1-choice);
//                	System.err.println("Grabbing random child.");
//                	logq -= Math.log(2); // TODO check ratio
//                }
//                //NodeRef oops
				//newParent =

//                logq += Math.log(possibleDestinations);

				tree.beginTreeEdit();

				// 4.1.1 if iP was root
				if (tree.isRoot(iP)) {
					// new root is CiP, but root cannot be a reassortment
					if (!tree.isBifurcation(CiP))
						throw new OperatorFailedException("root cannot be a reassortment");
					boolean doubly = tree.isBifurcationDoublyLinked(newParent);
					tree.doubleRemoveChild(iP, CiP);
					tree.doubleRemoveChild(newParent, newChild);

					//                   System.err.println("Down-slide disconnect:\n"+tree.toGraphString());
					if (tree.isBifurcation(newChild))
						tree.doubleAddChild(iP, newChild);
					else
						tree.singleAddChildWithOneParent(iP, newChild);
					if (!doubly)
						tree.doubleAddChild(newParent, iP);
					else {
						tree.singleAddChild(newParent, iP);
						tree.singleAddChildWithOneParent(newParent, newChild);
					}
					//tree.doubleAddChild(newParent, iP);
					//tree.swapHeightParameters(iP,CiP);  // TODO
					tree.setRoot(CiP);
//					System.err.println("DOWN: Creating new root!"); // TODO still not tested
				} else {
					boolean doubly = tree.isBifurcationDoublyLinked(newParent); // or is it CiP = newParent?
					tree.doubleRemoveChild(iP, CiP);
					tree.doubleRemoveChild(PiP, iP);
					tree.doubleRemoveChild(newParent, newChild);
					//               	System.err.println("Down-slide disconnect:\n"+tree.toGraphString());
					if (tree.isBifurcation(newChild))
						tree.doubleAddChild(iP, newChild);
					else
						tree.singleAddChildWithOneParent(iP, newChild);
					tree.doubleAddChild(PiP, CiP);
					if (!doubly)
						tree.doubleAddChild(newParent, iP); // Only works if iP is bifurcation
					else {
						tree.singleAddChild(newParent, iP);
						tree.singleAddChildWithOneParent(newParent, newChild);
					}
//					System.err.println("DOWN: no new root!");
					/*   	if( true ) {

																	 try {
																		 tree.endTreeEdit();
																	 } catch(MutableTree.InvalidTreeException ite) {
																		 throw new RuntimeException(ite.toString());
																	 }
																	 throw new OperatorFailedException("");
																 }  */
				}

				tree.setNodeHeight(iP, newHeight);
//				System.err.println("i         = " + ((Node) i).number);
//				System.err.println("iP        = " + ((Node) iP).number);
//				System.err.println("CiP       = " + ((Node) CiP).number);
//				System.err.println("newChild  = " + ((Node) newChild).number + " "
//						+ tree.getNodeHeight(newChild));
//				System.err.println("newParent = " + ((Node) newParent).number + " "
//						+ tree.getNodeHeight(newParent));
//				if (oops != null)
//					System.err.println("oops      = " + ((Node) oops).number + " "
//							+ tree.getNodeHeight(oops));
//				System.err.println("newHeight = " + newHeight);
				//              System.err.println("After slide down:\n"+tree.toGraphString());

                tree.endTreeEdit();
				try {
					tree.checkTreeIsValid();
				} catch (MutableTree.InvalidTreeException ite) {
					throw new RuntimeException(ite.toString());
				}

				//logq = -Math.log((double) possibleDestinations);
				logq += Math.log((double) possibleDestinations);


			} else {
				try {
					tree.setNodeHeight(iP, newHeight);
				} catch (Exception e) {
//					System.err.println("iP =" + ((Node) iP).number);
//					System.err.println("newHeight =" + newHeight);
					//System.exit(-1);
				}
				logq = 0.0;
			}
			//          System.err.println("After slide down:\n"+tree.toGraphString());
//			System.err.println("Sanity check in down-slide.");
//			sanityCheck();
		}
/*
        if (swapRates) {
            NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                double tmp = tree.getNodeRate(i);
                tree.setNodeRate(i, tree.getNodeRate(j));
                tree.setNodeRate(j, tmp);
            }

        }

        if (swapTraits) {
            NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                double tmp = tree.getNodeTrait(i);
                tree.setNodeTrait(i, tree.getNodeTrait(j));
                tree.setNodeTrait(j, tmp);
            }

        }
*/

		// todo fix
		// Reject all trees in which the root is doubly-
		if (tree.isBifurcationDoublyLinked(tree.getRoot()))
			throw new OperatorFailedException("invalid slide");

		// todo -- check all ARGTree.Roots
		if (!tree.validRoot())
			throw new OperatorFailedException("Roots are invalid");

		if (logq == Double.NEGATIVE_INFINITY)
			throw new OperatorFailedException("invalid slide");
//		System.err.println("Ending Subtree Slide Operation.");
		//System.err.println("logq = "+logq);
		//  logq = 0;


		if (scaledDirichletBranches) {
			if (oldTreeHeight != tree.getNodeHeight(tree.getRoot()))
				throw new OperatorFailedException("Temporarily disabled."); // TODO calculate Hastings ratio
		}


		return logq;
	}

	private double getDelta() {
		if (!gaussian) {
			return (MathUtils.nextDouble() * size) - (size / 2.0);
		} else {
			return MathUtils.nextGaussian() * size;
		}
	}

	private int getSlideableSubtrees(ARGModel tree, ArrayList<NodeRef> potentials) {
		int count = 0;
		for (int i = 0, n = tree.getNodeCount(); i < n; i++) {
			NodeRef node = tree.getNode(i);
			if (!tree.isRoot(node) && tree.isBifurcation(node)
					&& tree.isBifurcation(tree.getParent(node))) {
				if (potentials != null)
					potentials.add(node);
				count++;
			}
		}
		return count;
	}

	private int intersectingEdges(ARGModel tree, NodeRef node, NodeRef parent, double height, ArrayList<NodeRef[]> directChildren) {

		// if( tree.isBifurcation(node) ) {
		//ree.getParent(node); tree.getMinParentNodeHeight(node);
		//NodeRef parent = tree.getParent(node);
		if (tree.getNodeHeight(parent) < height) return 0;

		if (tree.getNodeHeight(node) < height) {
			if (directChildren != null) {
				NodeRef[] addition = new NodeRef[2];
				addition[0] = parent;
				addition[1] = node;
				directChildren.add(addition);
			}
			//directChildren.add(node);
			return 1;
		}

		int count = 0;
		//for (int i = 0; i < tree.getChildCount(node); i++) {
		count += intersectingEdges(tree, tree.getChild(node, 0), node, height, directChildren);
		if (tree.isBifurcation(node))
			count += intersectingEdges(tree, tree.getChild(node, 1), node, height, directChildren);

		return count;
	}
	// Handle reassortment nodes
//        // Only the parental edge which routines flows down is important
//        NodeRef parent0 = tree.getParent(node,0);
//        NodeRef parent1 = tree.getParent(node,1);
//        double thisHeight = tree.getNodeHeight(node);
//        double height0 = tree.getNodeHeight(parent0);
//        double height1 = tree.getNodeHeight(parent1);
//
//        if( (height0 < height) && (height1 < height) ) return 0;
//        int count = 0;
//        if( (thisHeight < height) && (height0 >= height) ) {
//        	if( directChildren !=null ) {
//        		NodeRef[] addition = new NodeRef[2];
//        		addition[0] = parent0;
//        		addition[1] = node;
//        		directChildren.add(addition);
//        	}
//        		//directChildren.add(node);
//        	count += 1;
//        }
//        if( (thisHeight < height) && (height1 >= height) ) {
//        	if( directChildren != null ) {
//         		NodeRef[] addition = new NodeRef[2];
//        		addition[0] = parent1;
//        		addition[1] = node;
//        		directChildren.add(addition);
//        	}
//        		//directChildren.add(node);
//        	count += 1;
//        }
//        if( count > 0 )
//        	return count;
//
//        return intersectingEdges(tree, tree.getChild(node,0), height, directChildren);
//    }

	/**
	 * @return the other child of the given parent.
	 */
	private NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

		if (tree.getChild(parent, 0) == child) {
			return tree.getChild(parent, 1);
		} else {
			return tree.getChild(parent, 0);
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
			boolean scaledDirichletBranches = false;

			CoercionMode mode = CoercionMode.parseMode(xo);

			if (xo.hasAttribute(SWAP_RATES)) {
				swapRates = xo.getBooleanAttribute(SWAP_RATES);
			}
			if (xo.hasAttribute(SWAP_TRAITS)) {
				swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
			}

			if (xo.hasAttribute(DIRICHLET_BRANCHES)) {
				scaledDirichletBranches = xo.getBooleanAttribute(DIRICHLET_BRANCHES);
			}

			ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
			int weight = xo.getIntegerAttribute("weight");
			double size = xo.getDoubleAttribute("size");
			boolean gaussian = xo.getBooleanAttribute("gaussian");
			return new ARGSubtreeSlideOperator(treeModel, weight, size, gaussian, swapRates,
					swapTraits, scaledDirichletBranches, mode);
		}

		public String getParserDescription() {
			return "An operator that slides a subtree.";
		}

		public Class getReturnType() {
			return SubtreeSlideOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				AttributeRule.newIntegerRule("weight"),
				AttributeRule.newDoubleRule("size"),
				AttributeRule.newBooleanRule("gaussian"),
				AttributeRule.newBooleanRule(SWAP_RATES, true),
				AttributeRule.newBooleanRule(SWAP_TRAITS, true),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
				new ElementRule(ARGModel.class)
		};
	};

}
