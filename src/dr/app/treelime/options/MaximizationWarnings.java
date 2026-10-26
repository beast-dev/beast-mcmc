/*
 * MaximizationWarnings.java
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

import dr.app.treelime.types.TreePriorType;

import java.util.ArrayList;
import java.util.List;

/**
 * Set-ups that generate valid XML but do not do what the user probably expects. The Maximization tab shows them,
 * Generate asks before writing, and the XML carries them as comments.
 */
public final class MaximizationWarnings {

    private MaximizationWarnings() {
    }

    public static List<String> get(TreeLimeOptions options) {
        List<String> warnings = new ArrayList<String>();
        MaximizationOptions maximization = options.maximizationOptions;

        // node heights selected for optimization but held fixed, e.g. by tip-date sampling
        for (Parameter parameter : options.selectParameters()) {
            if (maximization.classify(parameter, false) == MaximizationOptions.Treatment.NODE_HEIGHTS) {
                PartitionOptions owner = GradientRegistry.findOwner(options, parameter);
                if (owner instanceof PartitionTreeModel) {
                    GradientRegistry.NodeHeightPlan plan =
                            GradientRegistry.nodeHeightPlan(options, (PartitionTreeModel) owner);
                    if (!plan.isOptimized()) {
                        warnings.add("All node heights of " + treeName((PartitionTreeModel) owner) +
                                " are held fixed at the starting tree: " + plan.fixedReason + ".");
                    }
                }
            }
        }

        // calibrations are priors: maximum likelihood ignores them
        if (maximization.objective == MaximizationOptions.Objective.ML) {
            List<String> calibrations = new ArrayList<String>();
            for (Parameter parameter : options.selectParameters()) {
                if ((parameter.isNodeHeight || parameter.isStatistic) &&
                        options.treeModelOptions.isNodeCalibrated(parameter)) {
                    calibrations.add(parameter.getName());
                }
            }
            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                if (prior.getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION) {
                    calibrations.add("the calibrated Yule tree prior");
                }
            }
            if (!calibrations.isEmpty()) {
                warnings.add("Calibration priors are not part of the maximum-likelihood objective, so these are " +
                        "ignored: " + String.join(", ", calibrations) + ". Choose the MAP objective to use them.");
            }
        }
        return warnings;
    }

    private static String treeName(PartitionTreeModel tree) {
        return tree.getPrefix().isEmpty() ? "the tree" : "tree " + tree.getName();
    }
}
