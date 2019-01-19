/*
 * TreeParameterModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.HashMap;
import java.util.Map;

/**
 * This class maintains a parameter of length equal to the number of nodes in the tree.
 * It can optionally include the root node. If the root node is not included then this
 * class will track tree changes that change the root node number and swap the parameter
 * values so that the parameter values stay with their node when changes to the tree
 * topology occur.
 *
 * @author Alexei Drummond
 */
public class TreeParameterModel extends AbstractModel implements TreeTrait<Double>, TreeDoubleTraitProvider {

    protected final MutableTreeModel tree;

    // The tree parameter;
    private final Parameter parameter;

    // These are stored in parameters for check pointing reasons...

    // the index of the root node.
    private final Parameter rootNodeNumber;

    // mappings from node number to parameter index and back
    private final Parameter nodeNumberToParameterIndex;
    private final Parameter parameterIndexToNodeNumber;


    private boolean includeRoot = false;

    private Intent intent;


    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot tree if the parameter includes a value associated with the root node.
     */
    public TreeParameterModel(MutableTreeModel tree, Parameter parameter, boolean includeRoot) {
        this(tree, parameter, includeRoot, Intent.NODE);
    }

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot tree if the parameter includes a value associated with the root node.
     */
    public TreeParameterModel(MutableTreeModel tree, Parameter parameter, boolean includeRoot, Intent intent) {

        super("treeParameterModel");
        this.tree = tree;
        this.parameter = parameter;

        this.includeRoot = includeRoot;

        this.intent = intent;

        int dim = parameter.getDimension();
        int parameterSize = getParameterSize();
        if (dim != parameterSize) {
            if (dim < 2) {
                // if the parameter is of size 1 or hasn't had its dimension set then expand to fit
                parameter.setDimension(parameterSize);
            } else {
                throw new IllegalArgumentException("dimension of parameter does not match tree branch count");
            }
        }

        addModel(tree);
        addVariable(parameter);

        nodeNumberToParameterIndex = new Parameter.Default(parameter.getId() + ".nodeNumberToParameterIndex", tree.getNodeCount());
        parameterIndexToNodeNumber = new Parameter.Default(parameter.getId() + ".parameterIndexToNodeNumber", parameterSize);

        int k = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (includeRoot || !tree.isRoot(node)) {
                nodeNumberToParameterIndex.setParameterValue(i, k);
                parameterIndexToNodeNumber.setParameterValue(k, node.getNumber());
                k++;
            } else {
                nodeNumberToParameterIndex.setParameterValue(i, -1); // set the root index to an illegal value
            }
        }

        addVariable(nodeNumberToParameterIndex);
        addVariable(parameterIndexToNodeNumber);

        rootNodeNumber = new Parameter.Default(parameter.getId() + ".rootNodeNumber");
        rootNodeNumber.setParameterValue(0, tree.getRoot().getNumber());
        addVariable(rootNodeNumber);
    }

    public int getParameterSize() {
        return tree.getNodeCount() - (includeRoot ? 0 : 1);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            handleRootMove();
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == parameter) {
            // the underlying parameter has changed - fire an event for the associated node...
            int nodeNumber = getNodeNumberFromParameterIndex(index);

            assert (tree.getNode(nodeNumber).getNumber() == nodeNumber);

            fireModelChanged(variable, nodeNumber);
        }
    }

    protected void storeState() {
        // parameters store themselves
    }

    protected void restoreState() {
        // parameters restore themselves
    }

    protected void acceptState() {
    }

    public double getNodeDoubleValue(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    public double getNodeValue(Tree tree, NodeRef node) {

        assert (!tree.isRoot(node) || includeRoot) : "root node doesn't have a parameter value!";

        return parameter.getParameterValue(getParameterIndexFromNodeNumber(node.getNumber()));
    }

    public void setNodeValue(Tree tree, NodeRef node, double value) {

        assert (!tree.isRoot(node) && !includeRoot) : "root node doesn't have a parameter value!";

        parameter.setParameterValue(getParameterIndexFromNodeNumber(node.getNumber()), value);
    }

    public int getNodeNumberFromParameterIndex(int parameterIndex) {
        return parameterIndexToNodeNumber.getValue(parameterIndex).intValue();
    }

    public int getParameterIndexFromNodeNumber(int nodeNumber) {
        return nodeNumberToParameterIndex.getValue(nodeNumber).intValue();
    }

    private void handleRootMove() {

        if (!includeRoot) {

            final int oldRootNodeNumber = rootNodeNumber.getValue(0).intValue();
            final int newRootNodeNumber = tree.getRoot().getNumber();

            if (oldRootNodeNumber != newRootNodeNumber) {
                parameterIndexToNodeNumber.setParameterValue(getParameterIndexFromNodeNumber(newRootNodeNumber), oldRootNodeNumber);

                double oldRootParameterIndex = nodeNumberToParameterIndex.getParameterValue(oldRootNodeNumber);

                assert(oldRootParameterIndex < 0); // should be -1

                nodeNumberToParameterIndex.setParameterValue(oldRootNodeNumber, nodeNumberToParameterIndex.getParameterValue(newRootNodeNumber));
                nodeNumberToParameterIndex.setParameterValue(newRootNodeNumber, oldRootParameterIndex);
                
                rootNodeNumber.setParameterValue(0, newRootNodeNumber);
            }
        }
    }

    /**
     * @return the tree model that this parameter is synchronized with
     */
    public MutableTreeModel getTreeModel() {
        return tree;
    }

    public String[] getNodeAttributeLabel() {
        return new String[]{};
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        return new String[]{};
    }

    public String getTraitName() {
        return parameter.getId();
    }

    public Intent getIntent() {
        return intent;
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public boolean getLoggable() {
        return true;
    }

    public Double getTrait(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    public String getTraitString(Tree tree, NodeRef node) {
        return Double.toString(getNodeValue(tree, node));
    }
}
