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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.Arrays;

public class Weights extends AbstractModel implements RandomField.WeightProvider {

    private final TreeModel tree;
    private BigFastTreeIntervals intervals;

    private boolean indicesKnown;

    private int[] indices;

    public Weights(TreeModel tree) {
        super("just work");
        this.tree = tree;
        intervals = new BigFastTreeIntervals((TreeModel) tree);
        addModel(intervals);

        indicesKnown = false;
    }

    public double weight(int index1, int index2) {

        if (!indicesKnown) {
            indices = getIndex();
            indicesKnown = true;
        }

        if (Math.abs(index1 - index2) != 1) {
            return 0.0;
        } else {
            index1 = indices[index1];
            index2 = indices[index2];


            if (index2 < index1) {
                int temp = index2;
                index2 = index1;
                index1 = temp;
            }

            return 2/(intervals.getInterval(index1) + intervals.getInterval(index2));
        }

    }


    @Override
    public int getDimension() {

        int dim = 0;

        for (int i = 1; i < intervals.getIntervalCount(); ++i) {
            if (intervals.getIntervalTime(i) != intervals.getIntervalTime(i-1)) {
                dim += 1;
            }
        }

        return dim;

    }


    private int[] getIndex() {
        int dim = getDimension();
        int[] index = new int[dim + 1];
        int j = 0;
        for (int i = 1; i < intervals.getIntervalCount(); ++i) {
            if (intervals.getIntervalTime(i) != intervals.getIntervalTime(i-1)) {
                index[j] = i - 1;
                j += 1;
            }
        }
        index[dim] = intervals.getIntervalCount() - 1;
        return index;
    }


//    @Override
//    public void addModelListener(ModelListener listener) {
//
//    }
//
//    @Override
//    public void removeModelListener(ModelListener listener) {
//
//    }











    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == intervals) {
            indicesKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        indicesKnown = false;
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
