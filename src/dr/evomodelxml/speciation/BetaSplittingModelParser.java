/*
 * BetaSplittingModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.BetaSplittingModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a SpeciationModel. Recognises YuleModel.
 */
public class BetaSplittingModelParser extends AbstractXMLObjectParser {

    public static final String BETA_SPLITTING_MODEL = "betaSplittingModel";
    public static final String PHI = "phi";
    public static final String TREE = "branchingTree";


    public String getParserName() {
        return BETA_SPLITTING_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(PHI);
        Parameter phiParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TREE);
        Tree tree = (Tree) cxo.getChild(Tree.class);

        return new BetaSplittingModel(phiParameter, tree);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "The beta-splitting family of tree branching models (Aldous, 1996;2001).";
    }

    public Class getReturnType() {
        return BetaSplittingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PHI,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, "A parameter that ranges from -infinity (comb-tree) to +infinity (balanced tree)"),
            new ElementRule(TREE,
                    new XMLSyntaxRule[]{new ElementRule(Tree.class)})
    };
}
