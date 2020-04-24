/*
 * NodeHeightProxyParameter.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightProxyParameter extends Parameter.Proxy {

    private TreeModel tree;
    private TreeParameterModel indexHelper;

    public NodeHeightProxyParameter(String name,
                                    TreeModel tree,
                                    boolean includeRoot) {
        super(name, includeRoot ? tree.getInternalNodeCount() : tree.getInternalNodeCount() - 1);
        this.tree = tree;
        this.indexHelper = new TreeParameterModel(tree,
                new Parameter.Default(includeRoot ? tree.getInternalNodeCount() : tree.getInternalNodeCount() - 1, 0.0),
                includeRoot);
    }

    private int getNodeNumber(int index) {
        return indexHelper.getNodeNumberFromParameterIndex(index + tree.getExternalNodeCount());
    }

    @Override
    public double getParameterValue(int dim) {
        return tree.getNodeHeight(tree.getNode(getNodeNumber(dim)));
    }

    @Override
    public void setParameterValue(int dim, double value) {
        tree.setNodeHeight(tree.getNode(getNodeNumber(dim)), value);
        tree.pushTreeChangedEvent(tree.getNode(getNodeNumber(dim)));
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        tree.setNodeHeightQuietly(tree.getNode(getNodeNumber(dim)), value);
    }

    @Override
    public void fireParameterChangedEvent() {
        tree.pushTreeChangedEvent();
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    private static final String NODE_HEIGHT_PARAMETER = "nodeHeightProxyParameter";
    private static final String INCLUDE_ROOT = "rootNode";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            boolean includeRoot = xo.getBooleanAttribute(INCLUDE_ROOT);
            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            return new NodeHeightProxyParameter(NODE_HEIGHT_PARAMETER, tree, includeRoot);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(INCLUDE_ROOT),
                    new ElementRule(TreeModel.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return NodeHeightProxyParameter.class;
        }

        @Override
        public String getParserName() {
            return NODE_HEIGHT_PARAMETER;
        }
    };
}
