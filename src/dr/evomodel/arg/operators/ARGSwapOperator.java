/*
 * ARGSwapOperator.java
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

package dr.evomodel.arg.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGModel.Node;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.logging.Logger;

/**
 * This method moves the arg model around.  Use of both the
 * reassortment and bifurcation modes, as well as the event operator,
 * satisfies irreducibility.
 *
 * @author ebloomqu
 * @author Marc A. Suchard
 */
public class ARGSwapOperator extends SimpleMCMCOperator {

	public static final String ARG_SWAP_OPERATOR = "argSwapOperator";
	public static final String SWAP_TYPE = "type";
	public static final String BIFURCATION_SWAP = "bifurcationSwap";
	public static final String REASSORTMENT_SWAP = "reassortmentSwap";
	public static final String DUAL_SWAP = "dualSwap";
	public static final String FULL_SWAP = "fullSwap";
	public static final String NARROW_SWAP = "narrowSwap";

	private ARGModel arg;
	private String mode;

	public ARGSwapOperator(ARGModel arg, String mode, int weight) {
		this.arg = arg;
		this.mode = mode;

		setWeight(weight);
	}


	public double doOperation() {

		if (mode.equals(NARROW_SWAP)) {
			return narrowSwap();
		}


		if ((mode.equals(REASSORTMENT_SWAP) || mode.equals(DUAL_SWAP)) &&
				arg.getReassortmentNodeCount() == 0) {
			return 0.0;
		}

		ArrayList<NodeRef> bifurcationNodes = new ArrayList<NodeRef>(arg.getNodeCount());
		ArrayList<NodeRef> reassortmentNodes = new ArrayList<NodeRef>(arg.getNodeCount());

		setupBifurcationNodes(bifurcationNodes);
		setupReassortmentNodes(reassortmentNodes);

		if (mode.equals(BIFURCATION_SWAP)) {
			return bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
		} else if (mode.equals(REASSORTMENT_SWAP)) {
			return reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
		} else if (mode.equals(DUAL_SWAP)) {
			reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
			return bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
		}

		bifurcationNodes.addAll(reassortmentNodes);

		Collections.sort(bifurcationNodes, NodeSorter);

		for (NodeRef x : bifurcationNodes) {
			if (arg.isBifurcation(x))
				bifurcationSwap(x);
			else
				reassortmentSwap(x);
		}

		return 0;
	}

	private double narrowSwap() {
		ArrayList<NarrowSwap> possibleSwaps = new ArrayList<NarrowSwap>(arg.getNodeCount());
		findAllNarrowSwaps(possibleSwaps);
		int possibleSwapsBefore = possibleSwaps.size();

		if (possibleSwapsBefore == 0)
			return 0;

		doNarrowSwap(possibleSwaps.get(MathUtils.nextInt(possibleSwaps.size())));

		possibleSwaps.clear();
		findAllNarrowSwaps(possibleSwaps);

		return Math.log((double) possibleSwapsBefore / possibleSwaps.size());
	}

	public int findAllNarrowSwaps(ArrayList<NarrowSwap> moves) {

		for (int i = 0, n = arg.getInternalNodeCount(); i < n; i++) {
			Node x = (Node) arg.getInternalNode(i);
			if (x.bifurcation && !x.isRoot() && x.leftParent.bifurcation) {
				NarrowSwap a = new NarrowSwap(x.leftChild, x, x.leftParent);
				NarrowSwap b = new NarrowSwap(x.rightChild, x, x.leftParent);

				if (a.isValid())
					moves.add(a);
				if (b.isValid())
					moves.add(b);
			}
		}
		return moves.size();
	}

	private void doNarrowSwap(NarrowSwap swap) {
		arg.beginTreeEdit();

		String before = arg.toARGSummary();

		if (swap.c == swap.pb) {
			Node c = (Node) swap.c;
			Node p = (Node) swap.p;
			Node gp = (Node) swap.gp;

			if (c.leftParent == p) {
				c.leftParent = gp;
				c.rightParent = p;
			} else {
				c.leftParent = p;
				c.rightParent = gp;
			}
		} else if (arg.getChild(swap.p, 0) == arg.getChild(swap.p, 1)) {
			Node p = (Node) swap.p;
			Node c = (Node) swap.c;

			if (MathUtils.nextBoolean())
				p.leftChild = c.leftParent = null;
			else
				p.rightChild = c.rightParent = null;
			arg.removeChild(swap.gp, swap.pb);

			arg.singleAddChild(swap.gp, swap.c);
			arg.singleAddChild(swap.p, swap.pb);
		} else {
			arg.removeChild(swap.gp, swap.pb);
			arg.removeChild(swap.p, swap.c);
			arg.singleAddChild(swap.gp, swap.c);
			arg.singleAddChild(swap.p, swap.pb);
		}

		assert nodeCheck() : swap + " " + before + " " + arg.toARGSummary();

                    arg.pushTreeChangedEvent(swap.gp);
                    arg.pushTreeChangedEvent(swap.p);

        arg.endTreeEdit();
		try {
			arg.checkTreeIsValid();
		} catch (MutableTree.InvalidTreeException ite) {
			System.out.println(swap);
			System.out.println(before);
			System.err.println(ite.getMessage());
			System.exit(-1);
		} catch (NullPointerException e) {
			System.out.println(swap);
			System.out.println(before);
			System.err.println(e.getMessage());
			System.exit(-1);
		}
	}

	private class NarrowSwap {
		public NodeRef c;
		public NodeRef p;
		public NodeRef gp;
		public NodeRef pb;

		public NarrowSwap(NodeRef child, NodeRef parent, NodeRef gParent) {
			c = child;
			p = parent;
			gp = gParent;
			pb = arg.getOtherChild(gParent, parent);
		}

		public boolean isValid() {
//			if (arg.getNodeHeight(pb) < arg.getNodeHeight(p))
//				return true;
//			return false;
                        return (arg.getNodeHeight(pb) < arg.getNodeHeight(p));
		}

		public String toString() {
			return "Child: " + c.toString() +
					", Parent: " + p.toString() +
					", G-parent: " + gp.toString() +
					", P-brother: " + pb.toString();
		}

	}


	private double bifurcationSwap(NodeRef x) {
		Node startNode = (Node) x;

//		Node keepChild = startNode.leftChild;
		Node moveChild = startNode.rightChild;

		if (MathUtils.nextBoolean()) {
//			keepChild = moveChild;
			moveChild = startNode.leftChild;
		}

		ArrayList<NodeRef> possibleNodes = new ArrayList<NodeRef>(arg.getNodeCount());

		findNodesAtHeight(possibleNodes, startNode.getHeight());

		assert !possibleNodes.contains(startNode);
		assert possibleNodes.size() > 0;


		Node swapNode = (Node) possibleNodes.get(MathUtils.nextInt(possibleNodes.size()));
		Node swapNodeParent = swapNode.leftParent;


		arg.beginTreeEdit();

		String before = arg.toARGSummary();


		if (swapNode.bifurcation) {
			swapNodeParent = swapNode.leftParent;

			arg.singleRemoveChild(startNode, moveChild);

			if (swapNodeParent.bifurcation) {
				arg.singleRemoveChild(swapNodeParent, swapNode);
				arg.singleAddChild(swapNodeParent, moveChild);
			} else {
				arg.doubleRemoveChild(swapNodeParent, swapNode);
				arg.doubleAddChild(swapNodeParent, moveChild);
			}

			arg.singleAddChild(startNode, swapNode);

		} else {
			boolean leftSide = true;
			boolean[] sideOk = {swapNode.leftParent.getHeight() > startNode.getHeight(),
					swapNode.rightParent.getHeight() > startNode.getHeight()};

			if (sideOk[0] && sideOk[1]) {
				if (MathUtils.nextBoolean()) {
					swapNodeParent = swapNode.rightParent;
					leftSide = false;
				}
			} else if (sideOk[1]) {
				swapNodeParent = swapNode.rightParent;
				leftSide = false;
			}

			//Double linked parents
			if (swapNode.leftParent == swapNode.rightParent) {
				arg.singleRemoveChild(startNode, moveChild);

				if (leftSide) {
					swapNode.leftParent = null;
					swapNodeParent.leftChild = null;
				} else {
					swapNode.rightParent = null;
					swapNodeParent.rightChild = null;
				}

				arg.singleAddChild(startNode, swapNode);
				arg.singleAddChild(swapNodeParent, moveChild);
			} else if (swapNode.leftParent == startNode || swapNode.rightParent == startNode) {
				arg.singleRemoveChild(startNode, moveChild);

				if (swapNodeParent.bifurcation) {
					arg.singleRemoveChild(swapNodeParent, swapNode);
					arg.singleAddChild(swapNodeParent, moveChild);
				} else {
					arg.doubleRemoveChild(swapNodeParent, swapNode);
					arg.doubleAddChild(swapNodeParent, moveChild);
				}

				if (startNode.leftChild == null)
					startNode.leftChild = swapNode;
				else
					startNode.rightChild = swapNode;

				if (swapNode.leftParent == null)
					swapNode.leftParent = startNode;
				else
					swapNode.rightParent = startNode;

			} else {
				arg.singleRemoveChild(startNode, moveChild);

				if (swapNodeParent.bifurcation) {
					arg.singleRemoveChild(swapNodeParent, swapNode);
					arg.singleAddChild(swapNodeParent, moveChild);
				} else {
					arg.doubleRemoveChild(swapNodeParent, swapNode);
					arg.doubleAddChild(swapNodeParent, moveChild);
				}
				arg.singleAddChild(startNode, swapNode);
			}
		}

              arg.pushTreeChangedEvent(); // TODO Send only changed nodes
//              arg.pushTreeChangedEvent(startNode);
//              arg.pushTreeChangedEvent(swapNodeParent);

		assert nodeCheck();

        arg.endTreeEdit();
		try {
			arg.checkTreeIsValid();
		} catch (MutableTree.InvalidTreeException ite) {
			System.out.println(before);
			System.err.println(ite.getMessage());
			System.exit(-1);
		}

		return 0;
	}

	private double reassortmentSwap(NodeRef x) {
		Node startNode = (Node) x;
		Node startChild = startNode.leftChild;

		ArrayList<NodeRef> possibleNodes = new ArrayList<NodeRef>(arg.getNodeCount());

		findNodesAtHeight(possibleNodes, startNode.getHeight());

		assert !possibleNodes.contains(startNode);
		assert possibleNodes.size() > 0;

		Node swapNode = (Node) possibleNodes.get(MathUtils.nextInt(possibleNodes.size()));

		Node swapParent;

		arg.beginTreeEdit();

		if (swapNode.bifurcation) {
			swapParent = swapNode.leftParent;

			arg.doubleRemoveChild(startNode, startChild);

			if (swapParent.bifurcation)
				arg.singleRemoveChild(swapParent, swapNode);
			else
				arg.doubleRemoveChild(swapParent, swapNode);

			arg.doubleAddChild(startNode, swapNode);


			if (startChild.bifurcation) {
				startChild.leftParent = swapParent;
				startChild.rightParent = swapParent;
			} else {
				if (startChild.leftParent == null) {
					startChild.leftParent = swapParent;
				} else {
					startChild.rightParent = swapParent;
				}
			}
			if (!swapParent.bifurcation) {
				swapParent.leftChild = startChild;
				swapParent.rightChild = startChild;
			} else {
				if (swapParent.leftChild == null) {
					swapParent.leftChild = startChild;
				} else {
					swapParent.rightChild = startChild;
				}
			}

		} else {

			boolean leftSide = true;
			boolean[] sideOk = {swapNode.leftParent.getHeight() > startNode.getHeight(),
					swapNode.rightParent.getHeight() > startNode.getHeight()};

			swapParent = swapNode.leftParent;

			if (sideOk[0] && sideOk[1]) {
				if (MathUtils.nextBoolean()) {
					leftSide = false;
					swapParent = swapNode.rightParent;
				}
			} else if (sideOk[1]) {
				leftSide = false;
				swapParent = swapNode.rightParent;
			}

			if (swapNode.leftParent == swapNode.rightParent) {
				arg.doubleRemoveChild(startNode, startChild);

				if (leftSide) {
					swapParent.leftChild = swapNode.leftParent = null;
					swapParent.leftChild = startChild;
					swapNode.leftParent = startNode;
				} else {
					swapParent.rightChild = swapNode.rightParent = null;
					swapParent.rightChild = startChild;
					swapNode.rightParent = startNode;
				}

				startNode.leftChild = startNode.rightChild = swapNode;

				if (startChild.bifurcation) {
					startChild.leftParent = startChild.rightParent = swapParent;
				} else {
					if (startChild.leftParent == null)
						startChild.leftParent = swapParent;
					else
						startChild.rightParent = swapParent;
				}
			} else {
				arg.doubleRemoveChild(startNode, startChild);

				if (swapParent.bifurcation)
					arg.singleRemoveChild(swapParent, swapNode);
				else
					arg.doubleRemoveChild(swapParent, swapNode);

				startNode.leftChild = startNode.rightChild = swapNode;

				if (leftSide)
					swapNode.leftParent = startNode;
				else
					swapNode.rightParent = startNode;

				if (swapParent.bifurcation) {
					if (swapParent.leftChild == null)
						swapParent.leftChild = startChild;
					else
						swapParent.rightChild = startChild;
				} else {
					swapParent.leftChild = swapParent.rightChild = startChild;
				}

				if (startChild.bifurcation) {
					startChild.leftParent = startChild.rightParent = swapParent;
				} else {
					if (startChild.leftParent == null)
						startChild.leftParent = swapParent;
					else
						startChild.rightParent = swapParent;
				}

			}

		}

		arg.pushTreeChangedEvent();  // TODO Limit tree hit

        arg.endTreeEdit();
		try {
			arg.checkTreeIsValid();
		} catch (MutableTree.InvalidTreeException ite) {
			System.err.println(ite.getMessage());
			System.exit(-1);
		}


		return 0;
	}

	private void setupBifurcationNodes(ArrayList<NodeRef> list) {
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			NodeRef x = arg.getNode(i);
			if (arg.isInternal(x) && arg.isBifurcation(x) && !arg.isRoot(x)) {
				list.add(x);
			}
		}
	}

	private void setupReassortmentNodes(ArrayList<NodeRef> list) {
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			NodeRef x = arg.getNode(i);
			if (arg.isReassortment(x)) {
				list.add(x);
			}
		}
	}


	private void findNodesAtHeight(ArrayList<NodeRef> x, double height) {
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			Node test = (Node) arg.getNode(i);
			if (test.getHeight() < height) {
				if (test.bifurcation) {
					if (test.leftParent.getHeight() > height) {
						x.add(test);
					}
				} else {
					if (test.leftParent.getHeight() > height) {
						x.add(test);
					}
					if (test.rightParent.getHeight() > height) {
						x.add(test);
					}
				}
			}
		}
	}


	public String getOperatorName() {
		return mode;
	}

	public String getPerformanceSuggestion() {
		return "";
	}

	private Comparator<NodeRef> NodeSorter = new Comparator<NodeRef>() {

		public int compare(NodeRef o1, NodeRef o2) {
			double[] heights = {arg.getNodeHeight(o1), arg.getNodeHeight(o2)};

			if (heights[0] < heights[1]) {
				return -1;
			} else if (heights[0] > heights[1]) {
				return 1;
			}

			return 0;
		}
	};

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserDescription() {
			return "Swaps nodes on a tree";
		}

		public Class getReturnType() {
			return ARGSwapOperator.class;
		}

		private String[] validFormats = {BIFURCATION_SWAP, REASSORTMENT_SWAP,
				DUAL_SWAP, FULL_SWAP, NARROW_SWAP};

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule(WEIGHT),
				new StringAttributeRule(SWAP_TYPE, "The mode of the operator",
						validFormats, false),
				new ElementRule(ARGModel.class),

		};

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			int weight = xo.getIntegerAttribute(WEIGHT);

			String mode = xo.getStringAttribute(SWAP_TYPE);

			Logger.getLogger("dr.evomodel").info("Creating ARGSwapOperator: " + mode);

			ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
			return new ARGSwapOperator(arg, mode, weight);
		}

		public String getParserName() {
			return ARG_SWAP_OPERATOR;
		}

	};

	public boolean nodeCheck() {
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			Node x = (Node) arg.getNode(i);

			if (x.leftParent != x.rightParent &&
					x.leftChild != x.rightChild) {
				return false;
			}
			if (x.leftParent != null) {
				if (x.leftParent.leftChild.getNumber() != i &&
						x.leftParent.rightChild.getNumber() != i)
					return false;
			}
			if (x.rightParent != null) {
				if (x.rightParent.leftChild.getNumber() != i &&
						x.rightParent.rightChild.getNumber() != i)
					return false;
			}
			if (x.leftChild != null) {
				if (x.leftChild.leftParent.getNumber() != i &&
						x.leftChild.rightParent.getNumber() != i)
					return false;
			}
			if (x.rightChild != null) {
				if (x.rightChild.leftParent.getNumber() != i &&
						x.rightChild.rightParent.getNumber() != i)
					return false;
			}
		}

		return true;
	}

}
