/*
 * Weights.java
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

package dr.inference.distribution;


import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

public class GriddedWeights extends AbstractModel implements RandomField.WeightProvider {

    private final List<TreeModel> treeModelList;
    private final boolean rescaleByRootHeight;
    private final Parameter gridPoints;


    public GriddedWeights(List<TreeModel> treeModelList, Parameter gridPoints,
                          boolean rescaleByRootHeight) {
        super("Weight Provider with grid points");
        this.treeModelList = treeModelList;
        this.rescaleByRootHeight = rescaleByRootHeight;
        this.gridPoints = gridPoints;
//        addModel(tree);
//        addVariable(gridPoints);
    }

    public double weight(int index1, int index2) {
        // TODO check for the case where two grid points are equal? this should not happen in the first place; create interface
        if (Math.abs(index1 - index2) != 1) {
            return 0.0;
        } else {
            double higherEnd;
            double lowerEnd;
            int maxIndex = Math.max(index1, index2);
            if (maxIndex >= gridPoints.getDimension()) {
                TreeModel tree = treeModelList.get(0);
                if (treeModelList.size() == 1 && gridPoints.getParameterValue(maxIndex - 1) < tree.getNodeHeight(tree.getRoot())) {
                    higherEnd = tree.getNodeHeight(tree.getRoot());
                } else {
                    higherEnd = gridPoints.getParameterValue(maxIndex - 1) +
                            (gridPoints.getParameterValue(maxIndex - 1) - gridPoints.getParameterValue(maxIndex - 2)) / 2.0;
                }
            } else {
                higherEnd = gridPoints.getParameterValue(maxIndex);
            }
            if (maxIndex == 1) {
                lowerEnd = 0.0; // TODO what if the first grid point is actually 0? This should not be in the first place
            } else {
                lowerEnd = gridPoints.getParameterValue(maxIndex - 2);
            }
            return 2/(higherEnd - lowerEnd) * getFieldScalar();
        }

    }

    private double getFieldScalar() {
        if (rescaleByRootHeight) {
            if (treeModelList.size() > 1) {
                double sumRootsHeight = 0.0;
                for(TreeModel tree : treeModelList) {
                    sumRootsHeight += tree.getNodeHeight(tree.getRoot());
                }
                return sumRootsHeight / treeModelList.size();
            } else {
                return treeModelList.get(0).getNodeHeight(treeModelList.get(0).getRoot());
            }
        }
        return 1.0;
    }


    @Override
    public int getDimension() {
        return gridPoints.getDimension();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public boolean isUsed() {
        return false;
    }

    @Override
    public String getId() {
        return null;
    }

    @Override
    public void setId(String id) {

    }
}
