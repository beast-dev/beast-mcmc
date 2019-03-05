/*
 * NodeHeightTransformParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;


import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "nodeHeightTransform";
    private static final String NODEHEIGHT = "nodeHeights";
    private static final String RATIO = "ratios";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(NODEHEIGHT);
        Parameter nodeHeightParameter = (Parameter) cxo.getChild(Parameter.class);

        XMLObject dxo = xo.getChild(RATIO);
        Parameter ratioParameter = (Parameter) dxo.getChild(Parameter.class);

        if (ratioParameter.getDimension() == 1) {
            ratioParameter.setDimension(nodeHeightParameter.getDimension());
        }
        ratioParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ratioParameter.getDimension()));

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        NodeHeightTransform nodeHeightTransform = new NodeHeightTransform(nodeHeightParameter, ratioParameter, tree, branchRateModel);
        nodeHeightTransform.transform(nodeHeightParameter.getParameterValues(), 0, nodeHeightParameter.getDimension());
        return nodeHeightTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(RATIO, Parameter.class, "The ratio parameter"),
                new ElementRule(NODEHEIGHT, Parameter.class, "The nodeHeight parameter"),
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NodeHeightTransform.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
