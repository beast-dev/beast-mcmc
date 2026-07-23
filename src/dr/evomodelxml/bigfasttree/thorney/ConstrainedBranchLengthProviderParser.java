/*
 * ConstrainedBranchLengthProviderParser.java
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

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.IntegerDataType;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.MutationBranchMap;
import dr.xml.*;

public class ConstrainedBranchLengthProviderParser extends AbstractXMLObjectParser {

    public static final String MUTATION_BRANCH_MAP = "simpleMutationBranchMap";
    public static final String DATA_TREE = "dataTree";

    public String getParserName() {
        return MUTATION_BRANCH_MAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);
        boolean discrete = xo.getAttribute("discrete", true);
        double minBranchLength=xo.getAttribute("minBranchlength", 0.0);

        ConstrainedTreeModel constrainedTree = (ConstrainedTreeModel) xo.getChild(ConstrainedTreeModel.class);
        Tree dataTree = (Tree) xo.getElementFirstChild(DATA_TREE);
        DataType dataType = discrete? IntegerDataType.INSTANCE: ContinuousDataType.INSTANCE;

        return new ConstrainedTreeBranchLengthProvider(constrainedTree,dataTree,scale,minBranchLength,dataType);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a mapping between branches in a  constrained tree and the number of mutations on the same branch in a fixed tree ";
    }

    public Class getReturnType() {
        return MutationBranchMap.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(ConstrainedTreeModel.class),
            new ElementRule(DATA_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
            AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the branchlengths in the data tree such as sequence length. default is 1"),
            AttributeRule.newDoubleRule("minBranchlength",true,"minimum branch length to be provided for branches not in data tree. default is 0.0"),
            AttributeRule.newBooleanRule("discrete",true,"should branchlengths be rounded to a discrete number of mutations? default is true")
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
