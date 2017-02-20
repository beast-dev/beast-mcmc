/*
 * ARGExchangeOperator.java
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
import dr.evomodel.operators.ExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Implements branch exchange operations.
 * There is a NARROW and WIDE variety.
 * The narrow exchange is very similar to a rooted-tree
 * nearest-neighbour interchange but with the restriction
 * that node height must remain consistent.
 * <p/>
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */
public class ARGExchangeOperator extends SimpleMCMCOperator {

	public static final String NARROW_EXCHANGE = "argNarrowExchange";
	public static final String WIDE_EXCHANGE = "argWideExchange";

	public static final int NARROW = 0;
	public static final int WIDE = 1;

	private static final int MAX_TRIES = 10000;

	private int mode = NARROW;
	private ARGModel tree;

	public ARGExchangeOperator(int mode, ARGModel tree, int weight) {
		this.mode = mode;
		this.tree = tree;
		setWeight(weight);
	}

	public double doOperation() {

		double logHastings = 0.0;
		int tipCount = tree.getExternalNodeCount();

		if (mode == NARROW) {
			if (tree.getReassortmentNodeCount() < 2)
				logHastings = narrow();
			else
				return 0.0;
		} else {
			logHastings = wide();
		}


		if (tree.getExternalNodeCount() != tipCount) {
			throw new RuntimeException("Lost some tips in " + ((mode == NARROW) ? "NARROW mode." : "WIDE mode."));
		}

		assert !Double.isInfinite(logHastings) && !Double.isNaN(logHastings);

		return logHastings;
	}

	public int getAllValidNarrowMoves() {
		NodeRef iP = null, j = null, jP = null;
		ArrayList<NodeRef> nodes = new ArrayList<NodeRef>(tree.getNodeCount());
		ArrayList<NarrowMove> moves = new ArrayList<NarrowMove>(tree.getNodeCount());

		for (int k = 0, n = tree.getNodeCount(); k < n; k++) {
			NodeRef x = tree.getNode(k);
			if (!tree.isRoot(x) && !tree.isRoot(tree.getParent(x, 0))
					&& !tree.isRoot(tree.getParent(x, 1))) {
				nodes.add(x);
			}
		}
		NarrowMove a;
		for (NodeRef i : nodes) {
			for (int k = 0; k < 2; k++) {
				iP = tree.getParent(i, k);
				for (int m = 0; m < 2; m++) {
					jP = tree.getParent(iP, m);
					j = tree.getOtherChild(jP, iP);
					a = new NarrowMove(i, iP, j, jP);
					if (validMove(a) && !moves.contains(a))
						moves.add(a);
				}
			}
		}
		assert moves.size() > 0;
		return moves.size();
	}

	private boolean validMove(NarrowMove move) {
		if (move.j != move.iP && move.i != move.j &&
				(tree.getNodeHeight(move.j) < tree.getNodeHeight(move.iP)) &&
				(tree.getNodeHeight(move.i) < tree.getNodeHeight(move.jP))) {
			return true;

		}
		return false;
	}

	private class NarrowMove {
		public NodeRef i;
		public NodeRef j;
		public NodeRef iP;
		public NodeRef jP;

		public NarrowMove(NodeRef i, NodeRef iP, NodeRef j, NodeRef jP) {
			this.i = i;
			this.j = j;
			this.iP = iP;
			this.jP = jP;

		}

		public boolean equals(Object o) {
			if (!(o instanceof NarrowMove)) {
				return false;
			}

			NarrowMove move = (NarrowMove) o;

			if (this.i == move.i && this.j == move.j &&
					this.iP == move.iP && this.jP == move.jP) {
				return true;
			}
			if (this.i == move.j && this.j == move.i &&
					this.iP == move.jP && this.jP == move.iP) {
				return true;
			}

			return false;
		}

		public String toString() {
			return "(" + i.toString() + ", " + iP.toString() +
					", " + jP.toString() + ", " + j.toString() + ")";
		}

	}

	/**
	 * WARNING: Assumes strictly bifurcating tree.
	 */
	public double narrow() {

		NodeRef i = null, iP = null, j = null, jP = null;
		int tries = 0;

		//Echoose

		int beforeMoves = getAllValidNarrowMoves();

		while (tries < MAX_TRIES) {
			i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			while (tree.getRoot() == i || tree.getParent(i, 0) == tree.getRoot() || tree.getParent(i, 1) == tree.getRoot()) {
				i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			}

			iP = tree.getParent(i, 0);
			if (tree.isReassortment(i) && MathUtils.nextBoolean())
				iP = tree.getParent(i, 1);

			jP = tree.getParent(iP, 0);
			if (tree.isReassortment(iP) && MathUtils.nextBoolean())
				jP = tree.getParent(iP, 1);

			j = tree.getChild(jP, 0);
			if (j == iP) {
				j = tree.getChild(jP, 1);
			}

			if (j != iP && i != j &&             // can still occur if i is child of doubly-linked reassortment
					(tree.getNodeHeight(j) < tree.getNodeHeight(iP)) && (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
				// todo fix height check for cases where i and j get switched
				break;
			}
			tries += 1;
		}
		//System.out.println("tries = " + tries);

		//Eupdate
		if (tries < MAX_TRIES) {
			try {
				eupdateARG(i, j, iP, jP);
			} catch (ARGOperatorFailedException e) {
				return Double.NEGATIVE_INFINITY;
			}

			tree.pushTreeChangedEvent(iP);
			tree.pushTreeChangedEvent(jP);
		} else {
			//throw new ARGOperatorFailedException("Couldn't find valid narrow move on this tree!!");
			return Double.NEGATIVE_INFINITY;
		}

		return Math.log((double) beforeMoves / getAllValidNarrowMoves());
	}

	/**
	 * WARNING: Assumes strictly bifurcating tree.
	 */
	public double wide() {

		NodeRef i = null, iP = null, j = null, jP = null;
		int tries = 0;

		//Echoose

		while (tries < MAX_TRIES) {
			i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			while (tree.getRoot() == i) {
				i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			}

			j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			while (j == i || j == tree.getRoot()) {
				j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
			}

			iP = tree.getParent(i);
			jP = tree.getParent(j);

			if ((iP != jP) && (i != jP) && (j != iP) &&
					(tree.getNodeHeight(j) < tree.getNodeHeight(iP)) &&
					(tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
				break;
			}
			tries += 1;
		}
		//System.out.println("tries = " + tries);

		//Eupdate
		if (tries < MAX_TRIES) {
			try {
				eupdateARG(i, j, iP, jP);
			} catch (ARGOperatorFailedException e) {
				return Double.NEGATIVE_INFINITY;
			}
		} else {
//			throw new ARGOperatorFailedException("Couldn't find valid wide move on this tree!");
			return Double.NEGATIVE_INFINITY;
		}
		return 0.0;
	}

	public int getMode() {
		return mode;
	}

	public String getOperatorName() {
		return ((mode == NARROW) ? "Narrow" : "Wide") + " Exchange";
	}

	private void eupdateARG(NodeRef i, NodeRef j, NodeRef iP, NodeRef jP) throws ARGOperatorFailedException {

		// There are three different cases:
		// 1) neither i nor j are reassortments, 2) either i or j are reassortments, 3) both i and j are reassortments

		tree.beginTreeEdit();

		boolean iBifurcation = tree.isBifurcation(i);
		boolean jBifurcation = tree.isBifurcation(j);

		if (iBifurcation && jBifurcation) {
			tree.removeChild(iP, i);
			tree.removeChild(jP, j);
			tree.addChild(jP, i);
			tree.addChild(iP, j);
		} else if (!iBifurcation && !jBifurcation) {
//            tree.singleRemoveChild(iP, i);
//            tree.singleRemoveChild(jP, j);
//            tree.singleAddChild(jP, i);
//            tree.singleAddChild(iP, j);
		} else {
			if (jBifurcation) {
				// one reassortment; force i to be bifurcation and j to be reassortment
				NodeRef t = i;
				NodeRef tP = iP;
				i = j;
				iP = jP;
				j = t;
				jP = tP;
			}
//            System.err.println(tree.toGraphString());
//            ARGModel.Node iNode = (ARGModel.Node) i;
//            ARGModel.Node jNode = (ARGModel.Node) j;
//            System.err.println("i = "+iNode.number+" : j = "+jNode.number);
			/* tree.removeChild(iP, i);
					   tree.singleRemoveChild(jP, j);
					   tree.addChild(jP, i);
					   tree.singleAddChild(iP, j);*/
		}

        tree.endTreeEdit();
		try {
            tree.checkTreeIsValid();
		} catch (MutableTree.InvalidTreeException ite) {
			throw new ARGOperatorFailedException(ite.toString());
		}
	}

	public double getMinimumAcceptanceLevel() {
		if (mode == NARROW) return 0.05;
		else return 0.01;
	}

	public double getMinimumGoodAcceptanceLevel() {
		if (mode == NARROW) return 0.05;
		else return 0.01;
	}

	public String getPerformanceSuggestion() {
		if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
			return "";
		} else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
			return "";
		} else {
			return "";
		}
	}

	public static XMLObjectParser NARROW_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return NARROW_EXCHANGE;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
			int weight = xo.getIntegerAttribute("weight");

			return new ARGExchangeOperator(NARROW, treeModel, weight);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents a narrow exchange operator. " +
					"This operator swaps a random subtree with its uncle.";
		}

		public Class getReturnType() {
			return ExchangeOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule("weight"),
				new ElementRule(ARGModel.class)
		};

	};

	public static XMLObjectParser WIDE_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return WIDE_EXCHANGE;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
			int weight = xo.getIntegerAttribute("weight");

			return new ARGExchangeOperator(WIDE, treeModel, weight);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents a wide exchange operator. " +
					"This operator swaps two random subtrees.";
		}

		public Class getReturnType() {
			return ExchangeOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule("weight"),
				new ElementRule(ARGModel.class)
		};

	};
}