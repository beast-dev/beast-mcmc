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
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "nodeHeightTransform";
    private static final String NODEHEIGHT = "nodeHeights";
    private static final String RATIO = "ratios";
    private static final String COALESCENT_INTERVAL = "coalescentIntervals";
    private static final String REAL_LINE = "realLine";
    private static final String WITH_ROOT = "withRoot";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(NODEHEIGHT);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        Parameter ratioParameter = null;
        if (xo.hasChildNamed(RATIO)) {
            ratioParameter = (Parameter) xo.getChild(RATIO).getChild(Parameter.class);
        }

        if (ratioParameter != null) {
            if (ratioParameter.getDimension() == 1) {
                ratioParameter.setDimension(tree.getInternalNodeCount() - 1);
            }
            ratioParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ratioParameter.getDimension()));
        }

        Parameter coalescentIntervals = null;
        GMRFSkyrideLikelihood skyrideLikelihood = null;
        if (xo.hasChildNamed(COALESCENT_INTERVAL)) {
            cxo = xo.getChild(COALESCENT_INTERVAL);
            skyrideLikelihood = (GMRFSkyrideLikelihood) cxo.getChild(GMRFSkyrideLikelihood.class);
        }

        boolean withRoot = xo.getBooleanAttribute(WITH_ROOT);

        Transform nodeHeightTransform;
        if (ratioParameter != null) {
            NodeHeightTransform transform = new NodeHeightTransform(ratioParameter, tree, branchRateModel, withRoot);
            if (xo.getChild(RATIO).getAttribute(REAL_LINE, false)) {

                List<Transform> transforms = new ArrayList<Transform>();
                if (withRoot) {
                    transforms.add(new Transform.LogTransform());
                }
                for (int i = 0; i < ratioParameter.getDimension(); i++) {
                    transforms.add(new Transform.LogitTransform());
                }
                nodeHeightTransform = new Transform.ComposeMultivariable(new Transform.Array(transforms, null), transform);
            } else {
                nodeHeightTransform = transform;
            }
        } else {
            Parameter nodeHeightParameter = (Parameter) cxo.getChild(Parameter.class);
            nodeHeightTransform = new NodeHeightTransform(nodeHeightParameter, tree, skyrideLikelihood);
            coalescentIntervals = ((NodeHeightTransform) nodeHeightTransform).getParameter();
            cxo = xo.getChild(COALESCENT_INTERVAL);
            coalescentIntervals.setId(cxo.getId());
            cxo.setNativeObject(coalescentIntervals);
        }

        return nodeHeightTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(
                        new ElementRule(RATIO, Parameter.class, "The ratio parameter"),
                        new AndRule(
                                new ElementRule(COALESCENT_INTERVAL, GMRFSkyrideLikelihood.class,
                                        "Construct a proxy parameter for coalescent intervals from the Skyride likelihood."),
                                new ElementRule(NODEHEIGHT, Parameter.class, "The nodeHeight parameter")
                        )
                ),
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class),
                AttributeRule.newBooleanRule(REAL_LINE, true),
                AttributeRule.newBooleanRule(WITH_ROOT),
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Transform.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
