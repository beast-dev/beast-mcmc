/*
 * SetHeightsAction.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

package dr.app.tools.treeannotator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.*;

public class SetHeightsAction implements CladeAction {
    // the smallest a height range can be before we consider it a single value
    double HEIGHT_EPSILON = 3e-8; // one second in decimal years

    // the heights to set for the root of the tree
    private final List<Double> rootHeights;

    // the number of heights that need to be present before we consider the range and HPD
    private final int filterCount;

    SetHeightsAction(List<Double> rootHeights, int filterCount) {
        this.rootHeights = rootHeights;
        this.filterCount = filterCount;
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        assert tree instanceof MutableTree;
        BiClade biclade = (BiClade)clade;

        setCladeHeights(biclade, biclade.getHeightValues());

        if (clade.getBestLeft() == null) {
            return;
        }

//        setCladeHeights((BiClade)biclade.getBestLeft(), biclade.getLeftHeightValues());
//        setCladeHeights((BiClade)biclade.getBestRight(), biclade.getRightHeightValues());

        if (tree.isRoot(node)) {
            setCladeHeights(biclade, rootHeights);
        }
    }

    @Override
    public boolean expectAllClades() {
        return true;
    }

    public void setCladeHeights(BiClade clade, List<Double> heights) {
        if (clade == null || heights.isEmpty()) {
            return;
        }
        double[] values = heights.stream().mapToDouble(Double::doubleValue).toArray();
        Double[] range = getRange(values);

        clade.setMeanHeight(DiscreteStatistics.mean(values));
        clade.setMedianHeight(DiscreteStatistics.median(values));

        if (clade.getSize() > 1 && Math.abs(range[0] - range[1]) > HEIGHT_EPSILON && heights.size() >= filterCount) {
            clade.setHeightRange(range);
            clade.setHeightHPD(getHPDs(0.95, values));
        }
    }

    private static Double[] getRange(double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        return new Double[] {min, max};
    }

    private static Double[] getHPDs(double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * (double) values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        return new Double[] {lower, upper};
    }

}
