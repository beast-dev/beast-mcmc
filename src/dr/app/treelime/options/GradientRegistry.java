/*
 * GradientRegistry.java
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * How the gradient of the tree likelihoods with respect to a parameter is computed. The prior side of the gradient
 * is always taken numerically from the prior compound, which is cheap and covers every prior, tree prior and clock
 * prior; only the tree-likelihood side, the expensive part, has analytic forms. When in doubt the answer is NUMERIC,
 * which is slower but never incomplete.
 */
public final class GradientRegistry {

    private GradientRegistry() {
    }

    public enum LikelihoodTerm {
        NONE("none", "does not enter the tree likelihoods"),
        SUBSTITUTION("analytic", "<branchSubstitutionParameterGradient>"),
        SITE_RATE("analytic", "<siteRateModelGradient>"),
        BRANCH_RATE("analytic", "<branchRateGradient>"),
        NUMERIC("numeric", "finite differences of the tree likelihoods");

        private final String kind;
        private final String description;

        LikelihoodTerm(String kind, String description) {
            this.kind = kind;
            this.description = description;
        }

        public boolean isAnalytic() {
            return kind.equals("analytic");
        }

        public String getDescription() {
            return description;
        }

        public String toString() {
            return kind;
        }
    }

    private static final List<String> SUBSTITUTION_PARAMETERS = Arrays.asList("kappa", "kappa1", "kappa2");
    private static final List<String> SITE_RATE_PARAMETERS = Arrays.asList("alpha", "pInv");

    /**
     * The tree-likelihood term from the parameter's owner and name. The generator downgrades analytic terms to
     * numeric when some partition needs a tree likelihood without gradients.
     */
    public static LikelihoodTerm likelihoodTerm(TreeLimeOptions options, Parameter parameter) {
        PartitionOptions owner = findOwner(options, parameter);
        String name = parameter.getBaseName();

        if (owner instanceof PartitionTreePrior) {
            return LikelihoodTerm.NONE;
        }
        if (owner instanceof PartitionSubstitutionModel) {
            if (SUBSTITUTION_PARAMETERS.contains(name)) {
                return LikelihoodTerm.SUBSTITUTION;
            }
            if (SITE_RATE_PARAMETERS.contains(name)) {
                return LikelihoodTerm.SITE_RATE;
            }
        }
        if (owner instanceof PartitionClockModel) {
            ClockType clockType = ((PartitionClockModel) owner).getClockType();
            if ((clockType == ClockType.STRICT_CLOCK && name.equals("clock.rate")) ||
                    (clockType == ClockType.HMC_CLOCK && name.equals(ClockType.HMC_CLOCK_BRANCH_RATES))) {
                return LikelihoodTerm.BRANCH_RATE;
            }
        }
        return LikelihoodTerm.NUMERIC;
    }

    /**
     * @return the model that holds the parameter, or null; parameters made with a prior helper do not record it
     */
    public static PartitionOptions findOwner(TreeLimeOptions options, Parameter parameter) {
        if (parameter.getOptions() != null) {
            return parameter.getOptions();
        }
        List<PartitionOptions> models = new ArrayList<PartitionOptions>();
        models.addAll(options.getPartitionSubstitutionModels());
        models.addAll(options.getPartitionClockModels());
        models.addAll(options.getPartitionTreePriors());
        models.addAll(options.getPartitionTreeModels());
        for (PartitionOptions model : models) {
            if (model.parameters.get(parameter.getBaseName()) == parameter) {
                return model;
            }
        }
        return null;
    }

    /**
     * @return whether the tree likelihoods of these partitions depend on the parameter's owner
     */
    public static boolean usesOwner(TreeLimeOptions options, List<PartitionData> partitions, Parameter parameter) {
        PartitionOptions owner = findOwner(options, parameter);
        for (PartitionData partition : partitions) {
            if (partition.getPartitionSubstitutionModel() == owner || partition.getPartitionClockModel() == owner) {
                return true;
            }
        }
        return false;
    }
}
