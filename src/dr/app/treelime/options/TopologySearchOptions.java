/*
 * TopologySearchOptions.java
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

package dr.app.treelime.options;

import dr.app.treelime.types.ClockType;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * The Maximization tab's topology search (plan 7): off, or a hill-climb from the starting tree by rooted NNIs, or by
 * NNIs and clamped-height SPRs, or the Thorough preset (plan 7 M7, the most reliable at reaching the highest objective),
 * that accepts a change of topology only when the objective rises. Its round re-fits move the node heights and the branch
 * rates; everything else is re-fitted by the maximizer at the start, on a new best and at the end.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TopologySearchOptions implements Serializable {

    private static final long serialVersionUID = 1L;

    public enum Mode {
        // in the combo's order; files store the names
        OFF("Off: keep the starting tree's topology"),
        THOROUGH("Thorough: NNI, SPR and swaps from the NJ tree, then perturbations (best tree)"),
        NNI("NNI hill-climb from the starting tree"),
        NNI_SPR("NNI and SPR hill-climb from the starting tree"),
        // after the others: files store the names, and the mode is not in the GUI (plan 7 M7)
        STANDARD("Standard: NNI, SPR and root moves from the NJ tree (no perturbations)");

        private final String label;

        Mode(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }
    }

    public Mode mode = Mode.OFF;
    public long seed = 1;
    public String maxSearchTime = ""; // a duration with an optional unit, "" for none
    // perturbed hill-climbs in a row without a new best before stopping; 0: none, or the Thorough preset's 20 (0 in
    // older files)
    public int nStop = 0;
    // independent starts written by "Generate starts..." (separate BEAST runs); 0 in older files
    private int starts = 4;

    public int getStarts() {
        return starts >= 2 ? starts : 4;
    }

    public void setStarts(int starts) {
        this.starts = starts;
    }

    public boolean isOn() {
        return mode != null && mode != Mode.OFF;
    }

    /**
     * @return whether the search starts from the NJ tree of the data it is given (Standard and Thorough), not from
     * the starting tree; a searched bootstrap needs it, so that each replicate starts from its own data
     */
    public boolean startsFromNj() {
        return mode == Mode.STANDARD || mode == Mode.THOROUGH;
    }

    /**
     * @return why the search cannot run with these options, or null when it can
     */
    public static String getBlocker(TreeLimeOptions options) {
        MaximizationOptions maximization = options.maximizationOptions;
        if (maximization.objective == MaximizationOptions.Objective.PRIOR) {
            return "the prior-only objective does not depend on the data, so it cannot rank topologies";
        }
        if (!options.hasData()) {
            return "there are no data";
        }
        if (!maximization.isAnalytic()) {
            return "the search needs the analytic node-height gradients: choose analytic gradients";
        }
        if (options.getPartitionTreeModels().size() != 1) {
            return "the search handles one tree; these partitions have " + options.getPartitionTreeModels().size() +
                    " unlinked trees";
        }
        PartitionTreeModel tree = options.getPartitionTreeModels().get(0);
        if (getOptimizedHeights(options) != tree) {
            GradientRegistry.NodeHeightPlan plan = GradientRegistry.nodeHeightPlan(options, tree);
            return "the search needs the tree's node heights optimized" +
                    (plan.isOptimized() ? ": tick its nodeHeights row" : ", but they are " + plan.getDescription());
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            ClockType type = clock.getClockType();
            switch (type) {
                case STRICT_CLOCK:
                case HMC_CLOCK:
                case UNCORRELATED:
                    break;
                default:
                    return "the search does not support the " + type.toString().toLowerCase() +
                            ": its rates are tied to clades or to parent branches, which moves would break";
            }
            if (type != ClockType.STRICT_CLOCK && maximization.objective == MaximizationOptions.Objective.ML) {
                return "maximum likelihood with a relaxed clock does not identify the root or the node times " +
                        "(every rooting and dating gives the same branch lengths): use MAP or a strict clock";
            }
        }
        return null;
    }

    /**
     * @return the tree whose node heights the maximizer optimizes, or null for none (or several)
     */
    public static PartitionTreeModel getOptimizedHeights(TreeLimeOptions options) {
        MaximizationOptions maximization = options.maximizationOptions;
        List<PartitionTreeModel> trees = new ArrayList<PartitionTreeModel>();
        for (Parameter parameter : options.selectParameters()) {
            if (maximization.classify(parameter, false) == MaximizationOptions.Treatment.NODE_HEIGHTS) {
                PartitionOptions owner = GradientRegistry.findOwner(options, parameter);
                if (owner instanceof PartitionTreeModel &&
                        GradientRegistry.nodeHeightPlan(options, (PartitionTreeModel) owner).isOptimized()) {
                    trees.add((PartitionTreeModel) owner);
                }
            }
        }
        return trees.size() == 1 ? trees.get(0) : null;
    }
}
