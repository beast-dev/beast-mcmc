/*
 * DiscretizedBranchRates.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class LocalBranchRates extends AbstractModel implements BranchRateModel  {

    public static final String LOCAL_BRANCH_RATES = "localBranchRates";
    public static final String RATE_INDICATORS = "rateIndicator";
    public static final String RATES = "rates";

    // The rate categories of each branch

    // the index of the root node.
    private int rootNodeNumber;
    private int storedRootNodeNumber;

	public LocalBranchRates(TreeModel tree, Parameter rateIndicatorParameter, Parameter ratesParameter) {

        super(LOCAL_BRANCH_RATES);

        if (rateIndicatorParameter.getDimension() != tree.getNodeCount() -1 ) {
            throw new IllegalArgumentException("The rate category parameter must be of length nodeCount-1");
        }

        for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
            rateIndicatorParameter.setParameterValue(i, 0.0);
            ratesParameter.setParameterValue(i, 1.0);
        }

        addModel(tree);

        addParameter(rateIndicatorParameter);
        addParameter(ratesParameter);

        rootNodeNumber = tree.getRoot().getNumber();
        storedRootNodeNumber = rootNodeNumber;
	}

	public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        fireModelChanged();
    }

    protected void storeState() {
        storedRootNodeNumber = rootNodeNumber;
    }

    protected void restoreState() {
        rootNodeNumber = storedRootNodeNumber;
    }

    protected void acceptState() {}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return LOCAL_BRANCH_RATES; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel tree = (TreeModel)xo.getChild(TreeModel.class);

            Parameter rateIndicatorParameter = (Parameter)xo.getSocketChild(RATE_INDICATORS);
            Parameter ratesParameter = (Parameter)xo.getSocketChild(RATES);

			Logger.getLogger("dr.evomodel").info("Using local relaxed clock model.");

			return new LocalBranchRates(tree, rateIndicatorParameter, ratesParameter);
        }

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return
				"This element returns an discretized relaxed clock model." +
                "The branch rates are drawn from a discretized parametric distribution.";
		}

		public Class getReturnType() { return LocalBranchRates.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(TreeModel.class),
            new ElementRule(RATE_INDICATORS, Parameter.class, "The rate change indicators parameter", false),
            new ElementRule(RATES, Parameter.class, "The rate changes parameter", false)
        };
    };

    public double getBranchRate(Tree tree, NodeRef node) {

        if (tree.isRoot(node)) {
            return 1.0;
        } else {

            double rate;
            if (isRateChangeOnBranchAbove(tree, node)) {
                rate = tree.getNodeRate(node);
            } else {
                rate = getBranchRate(tree, tree.getParent(node));
            }
            return rate;
        }
    }

    public final boolean isRateChangeOnBranchAbove(Tree tree, NodeRef node) {
        return (int)Math.round(((TreeModel)tree).getNodeTrait(node)) == 1;
    }
}