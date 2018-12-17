/*
 * LocalBranchRatesParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.ArbitraryBranchRates.BranchRateTransform;
import dr.evomodel.branchratemodel.LocalBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class LocalBranchRatesParser extends AbstractXMLObjectParser {

    public static String SHRINKAGE_BRANCH_RATES = "localBranchRates";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        BranchRateTransform transform = (BranchRateTransform) xo.getChild(BranchRateTransform.class);

        if (transform == null) {
            transform = new BranchRateTransform.None();
        }

        Parameter multipllierParameter = (Parameter) xo.getChild(Parameter.class);

        final int numBranches = tree.getNodeCount() - 1;
        if (multipllierParameter.getDimension() != numBranches) {
            multipllierParameter.setDimension(numBranches);
        }

        return new LocalBranchRates(tree, multipllierParameter, transform);
    }

    @Override
    public String getParserDescription() {
        return "This element returns a local branch rate model. " +
                "The branch rates are drawn from shrinkage priors.";
    }

    @Override
    public Class getReturnType() {
        return LocalBranchRates.class;
    }

    @Override
    public String getParserName() {
        return SHRINKAGE_BRANCH_RATES;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(Parameter.class, "The multiplier parameter"),
            new ElementRule(ArbitraryBranchRates.BranchRateTransform.class, true),
    };
}
