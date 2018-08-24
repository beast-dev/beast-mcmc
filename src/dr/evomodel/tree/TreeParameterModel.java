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

    // the index of the root node.
    private final Parameter rootNodeNumber;

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
        int treeSize = getParameterSize();
        if (dim != treeSize) {
//            System.err.println("WARNING: setting dimension of parameter to match tree branch count ("
//                    + dim + " != " + treeSize + ")"); // http://code.google.com/p/beast-mcmc/issues/detail?id=385
            parameter.setDimension(treeSize);
        }

        addModel(tree);
        addVariable(parameter);

        rootNodeNumber = new Parameter.Default(parameter.getId() + ".rootNodeNumber");
        rootNodeNumber.setParameterValue(0, tree.getRoot().getNumber());
        addVariable(rootNodeNumber);
    }

    public int getParameterSize() {
        int treeSize = tree.getNodeCount();
        if (!includeRoot) {
            treeSize -= 1;
        }
        return treeSize;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            handleRootMove();
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        int nodeNumber = getNodeNumberFromParameterIndex(index);

        assert (tree.getNode(nodeNumber).getNumber() == nodeNumber);

        fireModelChanged(variable, nodeNumber);
    }

    protected void storeState() {
        //rootNodeNumber.storeParameterValues();
    }

    protected void restoreState() {
        //rootNodeNumber.restoreParameterValues();
    }

    protected void acceptState() {
    }

    public double getNodeDoubleValue(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    public double getNodeValue(Tree tree, NodeRef node) {

        assert (!tree.isRoot(node) || includeRoot) : "root node doesn't have a parameter value!";

        assert !includeRoot || tree.getRoot().getNumber() == rootNodeNumber.getValue(0).intValue() :
                "INTERNAL ERROR! node with number " + rootNodeNumber + " should be the root node.";

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);
        return parameter.getParameterValue(index);
    }

    public void setNodeValue(Tree tree, NodeRef node, double value) {

        assert (!tree.isRoot(node) && !includeRoot) : "root node doesn't have a parameter value!";

        assert tree.getRoot().getNumber() == rootNodeNumber.getValue(0).intValue() :
                "INTERNAL ERROR! node with number " + rootNodeNumber + " should be the root node.";

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);
        parameter.setParameterValue(index, value);
    }

    public int getNodeNumberFromParameterIndex(int parameterIndex) {
        if (!includeRoot && parameterIndex >= tree.getRoot().getNumber()) return parameterIndex + 1;
        return parameterIndex;
    }

    public int getParameterIndexFromNodeNumber(int nodeNumber) {
        if (!includeRoot && nodeNumber > tree.getRoot().getNumber()) return nodeNumber - 1;
        return nodeNumber;
    }

    private void handleRootMove() {

        if (!includeRoot) {

            final int oldRootNodeNumber = rootNodeNumber.getValue(0).intValue();
            final int newRootNodeNumber = tree.getRoot().getNumber();

            if (oldRootNodeNumber != newRootNodeNumber) {
                if (oldRootNodeNumber > newRootNodeNumber) {

                    final double oldValue = parameter.getParameterValue(newRootNodeNumber);

                    final int end = Math.min(parameter.getDimension() - 1, oldRootNodeNumber);
                    for (int i = newRootNodeNumber; i < end; i++) {
                        parameter.setParameterValue(i, parameter.getParameterValue(i + 1));
                    }

                    parameter.setParameterValue(end, oldValue);

                } else if (oldRootNodeNumber < newRootNodeNumber) {

                    final int end = Math.min(parameter.getDimension() - 1, newRootNodeNumber);

                    final double oldValue = parameter.getParameterValue(end);

                    for (int i = end; i > oldRootNodeNumber; i--) {
                        parameter.setParameterValue(i, parameter.getParameterValue(i - 1));
                    }

                    parameter.setParameterValue(oldRootNodeNumber, oldValue);
                }
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
