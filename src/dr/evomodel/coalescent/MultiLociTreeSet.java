/*
 * MultiLociTreeSet.java
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

package dr.evomodel.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;

import java.util.Arrays;
import java.util.List;

/**
 * @author Joseph Heled
 *         Date: 19/11/2007
 */
public interface MultiLociTreeSet {
    /**
     *
     * @return  Number of independent loci in set.
     */
    int nLoci();

    /**
     *
     * @param nt  index of tree to return
     * @return the nt's loci (i.e. tree)
     */
    Tree getTree(int nt);

    /**
     *
     * @param nt index of tree to return
     * @return  Coalecsent intervals for nt's tree
     */
    TreeIntervals getTreeIntervals(int nt);

    /**
     *
     * @param nt
     * @return Population factor of nt's tree
     */
    double getPopulationFactor(int nt);

    void storeTheState();

    void restoreTheState();

    public class Default implements MultiLociTreeSet, ModelListener {
        private final List<TreeModel> trees;
        private final List<Double> factors;
        private final boolean[] dirty;
        private final boolean[] gotDirty;

        private final TreeIntervals[] intervals;

        public Default(List<TreeModel> trees, List<Double> popFactors) {
            this.trees = trees;
            this.factors = popFactors;

            for(TreeModel t : trees ) {
                t.addModelListener(this);
            }

            final int nt = trees.size();
            dirty = new boolean[nt];
            gotDirty = new boolean[nt];
            intervals = new TreeIntervals[nt];
            Arrays.fill(dirty, true);
            Arrays.fill(intervals, null);
        }

        public int nLoci() {
            return trees.size();
        }

        public Tree getTree(int nt) {
            return trees.get(nt);
        }

        public TreeIntervals getTreeIntervals(int nt) {
            if( dirty[nt] ) {
                intervals[nt] = new TreeIntervals(trees.get(nt));
                intervals[nt].setMultifurcationLimit(0);
                dirty[nt] = false;
            }
            return intervals[nt];
        }

        public double getPopulationFactor(int nt) {
            return factors.get(nt);
        }

        public void storeTheState() {
            Arrays.fill(gotDirty, false);
        }

        public void restoreTheState() {
            for(int nt = 0; nt < gotDirty.length; ++nt) {
                if( gotDirty[nt] ) {
                    dirty[nt] = true;
                }
            }
        }

        public void modelChangedEvent(Model model, Object object, int index) {
            assert model instanceof TreeModel;
            final int k = trees.indexOf((TreeModel)model);
            dirty[k] = gotDirty[k] = true;
        }

        public void modelRestored(Model model) {}
    }
}
