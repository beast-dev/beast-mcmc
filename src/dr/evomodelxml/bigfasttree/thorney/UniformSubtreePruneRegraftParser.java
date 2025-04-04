/*
 * UniformSubtreePruneRegraftParser.java
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

package dr.evomodelxml.bigfasttree.thorney;

import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeOperator;
import dr.evomodel.bigfasttree.thorney.UniformSubtreePruneRegraft;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;


public class UniformSubtreePruneRegraftParser extends AbstractXMLObjectParser {
    public static final String UNIFORM_SUBTREE_PRUNE_REGRAFT = "uniformSubtreePruneRegraft";
    @Override
    public Object parseXMLObject(XMLObject  xo) throws XMLParseException {
        TreeModel tree =  (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);


       UniformSubtreePruneRegraft op = new UniformSubtreePruneRegraft(tree,weight);
       if(tree instanceof ConstrainedTreeModel){
           return ConstrainedTreeOperator.parse((ConstrainedTreeModel) tree, weight, op,xo);
       }else{
           return op;
       }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "An operator that prunes a regrafts within a clade.";
    }

    @Override
    public Class getReturnType() {
        return AbstractTreeOperator.class;
    }

    @Override
    public String getParserName() {
        return UNIFORM_SUBTREE_PRUNE_REGRAFT;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };
}
