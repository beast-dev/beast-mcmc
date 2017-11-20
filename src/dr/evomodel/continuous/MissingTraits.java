/*
 * MissingTraits.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.tree.TreeModel;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public interface MissingTraits {

    public void handleMissingTips();

    public boolean isCompletelyMissing(int index);

    public boolean isPartiallyMissing(int index);

    void computeWeightedAverage(double[] meanCache, int meanOffset0, double precision0,
                                int meanOffset1, double precision1, int meanThisOffset, int dim);

    abstract class Abstract implements MissingTraits {

        protected static final boolean DEBUG = false;

        Abstract(MutableTreeModel treeModel, List<Integer> missingIndices, int dim) {
            this.treeModel = treeModel;
            this.dim = dim;
            this.missingIndices = missingIndices;

            completelyMissing = new boolean[treeModel.getNodeCount()];
            Arrays.fill(completelyMissing, 0, treeModel.getExternalNodeCount(), false);
            Arrays.fill(completelyMissing, treeModel.getExternalNodeCount(), treeModel.getNodeCount(), true); // All internal and root nodes are missing
        }

        final protected MutableTreeModel treeModel;
        final protected int dim;
        final protected List<Integer> missingIndices;
        final protected boolean[] completelyMissing;
    }

    public class CompletelyMissing extends Abstract {

        CompletelyMissing(MutableTreeModel treeModel, List<Integer> missingIndices, int dim) {
            super(treeModel, missingIndices, dim);
        }

        public void handleMissingTips() {
            for (Integer i : missingIndices) {
                int whichTip = i / dim;
                Logger.getLogger("dr.evomodel").info(
                        "\tMarking taxon " + treeModel.getTaxonId(whichTip) + " as completely missing");
                completelyMissing[whichTip] = true;
            }
        }

        public boolean isCompletelyMissing(int index) {
            return completelyMissing[index];
        }

        public boolean isPartiallyMissing(int index) {
            return false;
        }

        public void computeWeightedAverage(double[] meanCache,
                 int meanOffset0, double precision0,
                 int meanOffset1, double precision1,
                 int meanThisOffset, int dim) {
            IntegratedMultivariateTraitLikelihood.computeWeightedAverage(
                    meanCache, meanOffset0, precision0,
                    meanCache, meanOffset1, precision1,
                    meanCache, meanThisOffset, dim);
        }
    }

    public class PartiallyMissing extends Abstract {

        PartiallyMissing(TreeModel treeModel, List<Integer> missingIndices, int dim) {
            super(treeModel, missingIndices, dim);
        }

        public void handleMissingTips() {
            throw new RuntimeException("Not yet implemented");
        }

        public boolean isCompletelyMissing(int index) {
            throw new RuntimeException("Not yet implemented");
        }

        public boolean isPartiallyMissing(int index) {
            throw new RuntimeException("Not yet implemented");
        }

        public void computeWeightedAverage(double[] meanCache, int meanOffset0, double precision0, int meanOffset1, double precision1, int meanThisOffset, int dim) {
            throw new RuntimeException("Not yet implemented");
        }
    }
}
