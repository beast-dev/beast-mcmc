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
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
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
    private static final String COALESCENT_INTERVAL = "coalescentIntervals";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(NODEHEIGHT);
        Parameter nodeHeightParameter = (Parameter) cxo.getChild(Parameter.class);

        Parameter ratioParameter = null;
        if (xo.hasChildNamed(RATIO)) {
            ratioParameter = (Parameter) xo.getChild(RATIO).getChild(Parameter.class);
        }

        if (ratioParameter != null) {
            if (ratioParameter.getDimension() == 1) {
                ratioParameter.setDimension(nodeHeightParameter.getDimension());
            }
            ratioParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ratioParameter.getDimension()));
        }

        Parameter.Default coalescentIntervals = null;
        GMRFSkyrideLikelihood skyrideLikelihood = null;
        if (xo.hasChildNamed(COALESCENT_INTERVAL)) {
            cxo = xo.getChild(COALESCENT_INTERVAL);
            coalescentIntervals = (Parameter.Default) cxo.getChild(Parameter.class);
            skyrideLikelihood = (GMRFSkyrideLikelihood) cxo.getChild(GMRFSkyrideLikelihood.class);
            double[] intervals = skyrideLikelihood.getCoalescentIntervals();
            if (coalescentIntervals.getDimension() != intervals.length) {
                throw new RuntimeException("Coalescent interval parameter should have dimension of " + skyrideLikelihood.getPopSizeParameter().getDimension() + "!\n");
            }
            //coalescentIntervals.substituteBuffer(intervals); // TODO This probably won't work due to Parameter double-buffering
        }

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        NodeHeightTransform nodeHeightTransform;
        if (ratioParameter != null) {
            nodeHeightTransform = new NodeHeightTransform(nodeHeightParameter, ratioParameter, tree, branchRateModel);
        } else {
            nodeHeightTransform = new NodeHeightTransform(nodeHeightParameter, coalescentIntervals, tree, skyrideLikelihood);
        }
        nodeHeightTransform.transform(nodeHeightParameter.getParameterValues(), 0, nodeHeightParameter.getDimension());
        return nodeHeightTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(new ElementRule(RATIO, Parameter.class, "The ratio parameter"),
                        new ElementRule(COALESCENT_INTERVAL, new XMLSyntaxRule[]{new ElementRule(Parameter.class), new ElementRule(GMRFSkyrideLikelihood.class)})
                        ),
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
