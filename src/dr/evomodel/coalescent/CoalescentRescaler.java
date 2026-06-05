/*
 * CoalescentRescaler.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.demographicmodel.ConstantPopulationModel;
import dr.evomodel.coalescent.demographicmodel.RescaleAwareDemographic;
import dr.evomodel.tree.TreeModel;

/**
 * Simulates a set of coalescent intervals given a demographic model.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class CoalescentRescaler {

    private TreeModel tree;

    private ConstantPopulationModel constantPopulationModel;

    private RescaleAwareDemographic rescaleAwareDemographic;

    private BigFastTreeIntervals bigFastTreeIntervals;

    public CoalescentRescaler(TreeModel tree, ConstantPopulationModel constantPopulationModel, RescaleAwareDemographic rescaleAwareDemographic) {
        this.tree = tree;
        this.constantPopulationModel = constantPopulationModel;
        this.rescaleAwareDemographic = rescaleAwareDemographic;
        this.bigFastTreeIntervals = new BigFastTreeIntervals("constant.population.size.tree", tree);
    }

    public TreeModel rescaleTree() {
        double[] intensityTimeProducts = new double[bigFastTreeIntervals.getIntervalCount()];
        double[] intervals = new double[bigFastTreeIntervals.getIntervalCount()];
        for (int i = 0; i < bigFastTreeIntervals.getIntervalCount(); i++) {
            intervals[i] = bigFastTreeIntervals.getInterval(i);
            intensityTimeProducts[i] = intervals[i] / ((ConstantPopulation) constantPopulationModel.getDemographicFunction()).getN0();
        }
        double[] scaledIntervals = rescaleAwareDemographic.rescaleInterval(intervals, intensityTimeProducts);

        double height = 0;
        for (int i = 0; i < scaledIntervals.length; i++) {
            height += scaledIntervals[i];
            int[] intervalNodeNumbers = bigFastTreeIntervals.getNodeNumbersForInterval(i);
            tree.setNodeHeightQuietly(tree.getNode(intervalNodeNumbers[1]), height);
        }
        tree.pushTreeChangedEvent();

        return tree;
    }

}
