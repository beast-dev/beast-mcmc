/*
 * CharacterTreeModelParser.java
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

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses a characterTreeModel: a character tree that may be multifurcating and whose topology may be sampled.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeModelParser extends AbstractXMLObjectParser {

    public static final String VARIABLE = "variable";
    public static final String NEWICK = "newick";
    public static final String KEEP = "keep";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        DataType dataType = DataTypeUtils.getDataType(xo);
        if (dataType == null) {
            dataType = Nucleotides.INSTANCE;
        }
        final boolean variable = xo.getAttribute(VARIABLE, false);
        final String newick;
        if (xo.hasAttribute(NEWICK)) {
            if (xo.getChild(Tree.class) != null) {
                throw new XMLParseException("characterTreeModel: give the newick attribute or a newick element, not both");
            }
            newick = xo.getStringAttribute(NEWICK);
        } else {
            final Tree tree = (Tree) xo.getChild(Tree.class);
            if (tree == null) {
                throw new XMLParseException("characterTreeModel: the character tree is required (a newick element or attribute)");
            }
            newick = TreeUtils.newick(tree);
        }
        final Parameter keep = xo.hasChildNamed(KEEP) ? (Parameter) xo.getElementFirstChild(KEEP) : null;
        final String id = xo.hasId() ? xo.getId() : CharacterTreeModel.CHARACTER_TREE_MODEL;
        final CharacterTreeModel model;
        try {
            model = new CharacterTreeModel(id, newick, dataType, variable, keep);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException("characterTreeModel '" + id + "': " + e.getMessage());
        }
        Logger.getLogger("dr.evomodel").info("Character tree model '" + id + "': " + model.getContractedNewick()
                + (keep != null ? " topology sampled (keep " + keep.getId() + ")," : variable ? " variable binary tree," : " fixed,")
                + " " + model.getContractedInternalNodeCount() + " internal nodes");
        return model;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(VARIABLE, true),
            AttributeRule.newStringRule(NEWICK, true),
            AttributeRule.newStringRule(DataType.DATA_TYPE, true),
            new ElementRule(Tree.class, true),
            new ElementRule(DataType.class, true),
            new ElementRule(KEEP, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

    @Override
    public String getParserDescription() {
        return "A character tree of a hierarchical substitution model that may be multifurcating; with keep its topology is "
                + "sampled over all multifurcating trees";
    }

    @Override
    public Class getReturnType() {
        return CharacterTreeModel.class;
    }

    @Override
    public String getParserName() {
        return CharacterTreeModel.CHARACTER_TREE_MODEL;
    }
}
