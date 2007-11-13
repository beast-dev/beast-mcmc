/*
 * RateExchangeOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.*;
import dr.xml.*;
import dr.math.MathUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Implements the RateExchange move.
 *
 * @author Alexei Drummond
 *
 * @version $Id: RateExchangeOperator.java,v 1.3 2005/01/06 14:46:36 rambaut Exp $
 */
public class RateExchangeOperator extends SimpleMCMCOperator {

	public static final String RATE_EXCHANGE = "rateExchange";
    public static final String SWAP_TRAITS = "swapTraits";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_AT_ROOT = "swapAtRoot";
    public static final String MOVE_HEIGHT = "moveHeight";

	private final TreeModel tree;
    private final boolean swapRates;
    private final boolean swapTraits;
    private final boolean swapAtRoot;
    private final boolean moveHeight;

	public RateExchangeOperator(TreeModel tree, int weight, boolean swapRates, boolean swapTraits, boolean swapAtRoot, boolean moveHeight) {
		this.tree = tree;
		setWeight(weight);
        this.swapRates = swapRates;
        this.swapTraits = swapTraits;
        this.swapAtRoot = swapAtRoot;
        this.moveHeight = moveHeight;
	}

    /**
     * Do a probablistic subtree slide move.
     * @return the log-transformed hastings ratio
     */
	public double doOperation() throws OperatorFailedException {

        NodeRef node0 = tree.getInternalNode(MathUtils.nextInt(tree.getInternalNodeCount()));
        NodeRef node1 = tree.getChild(node0, 0);
        NodeRef node2 = tree.getChild(node0, 1);

        if (swapRates) {
            if (swapAtRoot) {
                double[] rates = new double[] { tree.getNodeRate(node0), tree.getNodeRate(node1), tree.getNodeRate(node2) };

                int r1 = MathUtils.nextInt(3);
                tree.setNodeRate(node0, rates[r1]);
                // swap down the top trait
                rates[r1] = rates[2];

                int r2 = MathUtils.nextInt(2);
                tree.setNodeRate(node1, rates[r2]);
                // swap down the top trait
                rates[r2] = rates[1];

                tree.setNodeRate(node2, rates[0]);
            } else {
                // just swap the two child rates...
                double tmp = tree.getNodeRate(node1);
                tree.setNodeRate(node1, tree.getNodeRate(node2));
                tree.setNodeRate(node2, tmp);
            }
        }

        if (swapTraits) {
            if (swapAtRoot) {
                double[] traits = new double[] { tree.getNodeTrait(node0), tree.getNodeTrait(node1), tree.getNodeTrait(node2) };

                int r1 = MathUtils.nextInt(3);
                tree.setNodeTrait(node0, traits[r1]);
                // swap down the top trait
                traits[r1] = traits[2];

                int r2 = MathUtils.nextInt(2);
                tree.setNodeTrait(node1, traits[r2]);
                // swap down the top trait
                traits[r2] = traits[1];

                tree.setNodeTrait(node2, traits[0]);
            } else {
                // just swap the two child traits...
                double tmp = tree.getNodeTrait(node1);
                tree.setNodeTrait(node1, tree.getNodeTrait(node2));
                tree.setNodeTrait(node2, tmp);
            }
        }

        // If the node is not the root, do a uniform pick of its height
        if (!tree.isRoot(node0) && moveHeight) {
            double lower = tree.getNodeHeightLower(node0);
            double upper = tree.getNodeHeightUpper(node0);
            double newValue = (MathUtils.nextDouble() * (upper - lower)) + lower;
            tree.setNodeHeight(node0, newValue);
        }

        return 0.0;
    }

	public double getTargetAcceptanceProbability() { return 0.234; }


	public String getPerformanceSuggestion() {
        if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()){
            return "";
        } else {
            return "";
        }
	}

	public String getOperatorName() {
		return RATE_EXCHANGE;
	}

	public Element createOperatorElement(Document d) {
		Element e = d.createElement(RATE_EXCHANGE);
		return e;
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() { return RATE_EXCHANGE; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
			int weight = xo.getIntegerAttribute("weight");
            boolean swapRates = xo.getBooleanAttribute(SWAP_RATES);
            boolean swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
            boolean swapAtRoot = xo.getBooleanAttribute(SWAP_AT_ROOT);
            boolean moveHeight = xo.getBooleanAttribute(MOVE_HEIGHT);
			return new RateExchangeOperator(treeModel, weight, swapRates, swapTraits, swapAtRoot, moveHeight);
		}

		public String getParserDescription() {
			return "An operator that exchanges rates and traits on a tree.";
		}

		public Class getReturnType() { return SubtreeSlideOperator.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newIntegerRule("weight"),
            AttributeRule.newBooleanRule(SWAP_RATES),
            AttributeRule.newBooleanRule(SWAP_TRAITS),
            AttributeRule.newBooleanRule(SWAP_AT_ROOT),
            AttributeRule.newBooleanRule(MOVE_HEIGHT),
			new ElementRule(TreeModel.class)
		};
	};

}
