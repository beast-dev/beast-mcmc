/*
 * AbstractTreeParameterModel.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
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
 * @author Marc A. Suchard
 */
public abstract class AbstractTreeParameterModel<T> extends AbstractModel implements TreeTrait<T> {

    protected final TreeModel tree;

    // The tree parameter;
    protected final Parameter parameter;

    // the index of the root node.
    private int rootNodeNumber;
    private int storedRootNodeNumber;

    private boolean includeRoot = false;
    private boolean includeTips = true;

    private final Intent intent;

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot true if the parameter includes values associated with the root node.
     * @param includeTips true if the parameter includes values associated with the tip nodes.
     */
    public AbstractTreeParameterModel(TreeModel tree, Parameter parameter, boolean includeRoot, boolean includeTips) {
        this(tree, parameter, includeRoot, includeTips, Intent.BRANCH);
    }

    public AbstractTreeParameterModel(TreeModel tree, Parameter parameter, boolean includeRoot, boolean includeTips,
                                      Intent intent) {

        super("treeParameterModel");
        this.tree = tree;
        this.parameter = parameter;

        this.includeRoot = includeRoot;
        this.includeTips = includeTips;

        this.intent = intent;

        int dim = parameter.getDimension();
        int treeSize = getParameterSize();
        if (dim != treeSize) {
            parameter.setDimension(treeSize);
        }

        addModel(tree);
        addVariable(parameter);

        rootNodeNumber = tree.getRoot().getNumber();
        storedRootNodeNumber = rootNodeNumber;
    }

    public abstract int getParameterSize();

    protected Parameter getParameter() {
        return parameter;
    }

    protected boolean checkNode(Tree tree, NodeRef node) {
        assert (!tree.isRoot(node) && !doesIncludeRoot()) : "root node doesn't have a parameter value!";

        assert tree.getRoot().getNumber() == rootNodeNumber :
                "INTERNAL ERROR! node with number " + rootNodeNumber + " should be the root node.";

        assert (!tree.isExternal(node) && !includeTips) : "tip nodes do not have parameter values!";

        return true;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            handleRootMove();
        }
    }

    protected boolean doesIncludeRoot() {
        return includeRoot;
    }

    protected boolean doesIncludeTips() {
        return includeTips;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        int nodeNumber = getNodeNumberFromParameterIndex(index);

        assert (tree.getNode(nodeNumber).getNumber() == nodeNumber);

        fireModelChanged(variable, nodeNumber);
    }

    protected void storeState() {
        storedRootNodeNumber = rootNodeNumber;
    }

    protected void restoreState() {
        rootNodeNumber = storedRootNodeNumber;
    }

    protected void acceptState() {
    }

    public abstract T getNodeValue(Tree tree, NodeRef node);

    public abstract void setNodeValue(Tree tree, NodeRef node, T value);

    protected int getNodeNumberFromParameterIndex(int parameterIndex) {
        int number = parameterIndex;
        if (!includeRoot && parameterIndex >= tree.getRoot().getNumber()) {
            number += 1;
        }
        if (!includeTips) {
            number += tree.getExternalNodeCount();
        }
        return number;
    }

    protected int getParameterIndexFromNodeNumber(int nodeNumber) {
        int index = nodeNumber;
        if (!includeRoot && nodeNumber > tree.getRoot().getNumber()) {
            index -= 1;
        }

        if (!includeTips) {
            index -= tree.getExternalNodeCount();
        }
        return index;

    }

    private void handleRootMove() {

        if (!includeRoot) {

            final int newRootNodeNumber = tree.getRoot().getNumber();

            if (rootNodeNumber > newRootNodeNumber) {

                final double oldValue = parameter.getParameterValue(newRootNodeNumber);

                final int end = Math.min(parameter.getDimension() - 1, rootNodeNumber);
                for (int i = newRootNodeNumber; i < end; i++) {
                    parameter.setParameterValue(i, parameter.getParameterValue(i + 1));
                }

                parameter.setParameterValue(end, oldValue);

            } else if (rootNodeNumber < newRootNodeNumber) {

                final int end = Math.min(parameter.getDimension() - 1, newRootNodeNumber);

                final double oldValue = parameter.getParameterValue(end);

                for (int i = end; i > rootNodeNumber; i--) {
                    parameter.setParameterValue(i, parameter.getParameterValue(i - 1));
                }

                parameter.setParameterValue(rootNodeNumber, oldValue);
            }
            rootNodeNumber = newRootNodeNumber;
        }
    }

    /**
     * @return the tree model that this parameter is synchronized with
     */
    public TreeModel getTreeModel() {
        return tree;
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

    public abstract T getTrait(Tree tree, NodeRef node);

    public abstract String getTraitString(Tree tree, NodeRef node);

}
