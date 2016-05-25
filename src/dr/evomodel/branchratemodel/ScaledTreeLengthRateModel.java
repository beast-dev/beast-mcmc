/*
 * ScaledTreeLengthRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.ScaledTreeLengthRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Package: ScaledTreeLengthRateModel
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 18, 2008
 * Time: 3:58:43 PM
 */
public class ScaledTreeLengthRateModel extends AbstractBranchRateModel {

    private Parameter totalLength;
    protected Tree treeModel;
    private double storedRateFactor;
    private boolean currentFactorKnown;
    private double rateFactor;

    public ScaledTreeLengthRateModel(TreeModel treeModel, Parameter totalLength) {
        super(ScaledTreeLengthRateModelParser.MODEL_NAME);
        this.totalLength = totalLength;
        this.treeModel = treeModel;
        currentFactorKnown = false;
        addModel(treeModel);
        addVariable(totalLength);
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        assert(tree == treeModel);

        if (!currentFactorKnown) {
            updateCurrentLength();
        }
        return rateFactor;
    }

    public double getTotalLength() {
        return totalLength.getParameterValue(0);
    }

    protected void updateCurrentLength() {
        double currentLength = 0;
        NodeRef root = treeModel.getRoot();
        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
            NodeRef node = treeModel.getNode(i);
            if (node != root) {
                currentLength += treeModel.getBranchLength(node);
            }
        }
        rateFactor = totalLength.getParameterValue(0) / currentLength;
        currentFactorKnown = true;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            currentFactorKnown = false;
        }
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == totalLength) {
            currentFactorKnown = false;
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
        storedRateFactor = rateFactor;
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
        rateFactor = storedRateFactor;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

}
