/*
 * BranchSpecificCompoundParameterParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.model.BranchParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchParameterParser extends AbstractXMLObjectParser {

    public static final String BRANCH_SPECIFIC_COMPOUND_PARAMETER = "branchParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rootParameter = (Parameter) xo.getChild(Parameter.class);
        ArbitraryBranchRates branchRates = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);
        Tree tree = branchRates.getTree();

        ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);

        BranchParameter branchParameter = new BranchParameter(
                xo.getName(),
                tree,
                branchRateModel,
                rootParameter
                );
        return branchParameter;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ArbitraryBranchRates.class),
                new ElementRule(Parameter.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "A multidimensional parameter constructed for each branch of the tree.";
    }

    @Override
    public Class getReturnType() {
        return BranchParameter.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_SPECIFIC_COMPOUND_PARAMETER;
    }
}
