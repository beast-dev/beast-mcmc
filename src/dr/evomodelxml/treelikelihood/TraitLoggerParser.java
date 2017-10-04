/*
 * StratifiedTraitLoggerParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treelikelihood.utilities.TreeTraitLogger;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class TraitLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "traitLogger";
    public static final String TRAIT_NAME = "traitName";
    private static final String NODES = "nodes";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree treeModel = (Tree) xo.getChild(Tree.class);
        TreeTrait trait = parseTreeTrait(xo, false);
        TreeTraitLogger.NodeRestriction nodes = TreeTraitLogger.NodeRestriction.parse(
                xo.getAttribute(NODES, "all"));

        return new TreeTraitLogger(treeModel, new TreeTrait[] { trait }, nodes);
    }

    static TreeTrait parseTreeTrait(XMLObject xo, boolean wholeTreeOnly) throws XMLParseException {

        TreeTraitProvider traitProvider =
                (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

        String traitName = (String) xo.getAttribute(TRAIT_NAME);

        TreeTrait trait = traitProvider.getTreeTrait(traitName);

        if (trait == null || (wholeTreeOnly && trait.getIntent() != TreeTrait.Intent.WHOLE_TREE)) {

            StringBuilder sb = new StringBuilder("Unable to find ");

            if (wholeTreeOnly) {
                sb.append("whole ");
            }

            String id = xo.hasId() ? xo.getId() : "null";

            sb.append("tree trait '").append(traitName).append("' in '").append(id).append("\n");

            sb.append("\tPossible traits:");
            for (TreeTrait existingTrait : traitProvider.getTreeTraits()) {
                if (!wholeTreeOnly || existingTrait.getIntent() == TreeTrait.Intent.WHOLE_TREE) {
                    sb.append(" ").append(existingTrait.getTraitName());
                }
            }
            sb.append("\n");

            throw new XMLParseException(sb.toString());
        }

        return trait;
    }

  private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
          AttributeRule.newStringRule(TRAIT_NAME),
          AttributeRule.newStringRule(NODES, true),
          new ElementRule(Tree.class),
          new ElementRule(TreeTraitProvider.class),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "A parser to log tree traits";
    }

    public Class getReturnType() {
        return TreeTraitLogger.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
