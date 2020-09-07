/*
 * BeagleTreeLikelihoodParser.java
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

package dr.evomodelxml.bigfasttree;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.ApproximatePoissonTreeLikelihood;
import dr.evomodel.bigfasttree.BranchLengthProvider;
import dr.evomodel.bigfasttree.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.bigfasttree.RzhetskyNeiBranchLengthProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.matrix.Matrix;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @author JT McCrone
 * @version $Id$
 */
public class ApproximatePoissonTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "approximatePoissonTreeLikelihood";
    public static final String DATA = "data";


    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int sequenceLength = xo.getIntegerAttribute("sequenceLength");

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        BranchLengthProvider branchLengthProvider;
        if(xo.getElementFirstChild(DATA) instanceof  Tree){
            Tree dataTree = (Tree) xo.getElementFirstChild(DATA);
            branchLengthProvider = new ConstrainedTreeBranchLengthProvider(dataTree,treeModel);
        }else{
            DistanceMatrix dataMatrix = (DistanceMatrix)xo.getElementFirstChild(DATA);
            branchLengthProvider = new RzhetskyNeiBranchLengthProvider(dataMatrix,treeModel);
        }

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        return new ApproximatePoissonTreeLikelihood(TREE_LIKELIHOOD, sequenceLength, treeModel, branchRateModel,branchLengthProvider);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule("sequenceLength", false),
            new OrRule(
                    new ElementRule(DATA, new XMLSyntaxRule[]{
                            new ElementRule(Tree.class)
                    }),
                    new ElementRule(DATA, new XMLSyntaxRule[]{
                            new ElementRule(DistanceMatrix.class)
                    })
            ),
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
