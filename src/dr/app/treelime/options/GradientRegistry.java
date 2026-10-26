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

import dr.app.treelime.components.tipdatesampling.TipDateSamplingComponentOptions;
import dr.app.treelime.types.ClockType;
import dr.app.treelime.types.PriorType;
import dr.app.treelime.types.TipDateSamplingType;
import dr.evolution.util.Taxa;

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
                    (clockType == ClockType.HMC_CLOCK && name.equals(ClockType.HMC_CLOCK_BRANCH_RATES)) ||
                    localClockIndex(options, (PartitionClockModel) owner, parameter) >= 0) {
                return LikelihoodTerm.BRANCH_RATE;
            }
        }
        return LikelihoodTerm.NUMERIC;
    }

    /**
     * The fixed local clock's branch-rate gradient is taken with respect to all its rates at once, stacked as the
     * clock writes them: the global rate, then one rate per monophyletic taxon set, in taxon-set order.
     *
     * @return the parameter's position in that stack, or -1 when it is not a rate of a fixed local clock
     */
    public static int localClockIndex(TreeLimeOptions options, PartitionClockModel clock, Parameter parameter) {
        if (clock.getClockType() != ClockType.FIXED_LOCAL_CLOCK) {
            return -1;
        }
        if (parameter == clock.getClockRateParameter()) {
            return 0;
        }
        int index = 1;
        for (Taxa taxonSet : options.taxonSets) {
            if (options.taxonSetsMono.get(taxonSet)) {
                if (clock.hasParameter(taxonSet.getId() + ".rate") &&
                        clock.getParameter(taxonSet.getId() + ".rate") == parameter) {
                    return index;
                }
                index++;
            }
        }
        return -1;
    }

    /**
     * @return the number of rates of a fixed local clock: the global rate and one per monophyletic taxon set
     */
    public static int localClockDimension(TreeLimeOptions options) {
        int dimension = 1;
        for (Taxa taxonSet : options.taxonSets) {
            if (options.taxonSetsMono.get(taxonSet)) {
                dimension++;
            }
        }
        return dimension;
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

    /**
     * How the internal node heights of one tree are optimized: held fixed (with the reason), or with the tree
     * likelihoods' analytic node-height gradient and, on the prior side, either the tree prior's analytic gradient
     * plus any CTMC rate reference priors on the tree, or finite differences of the prior compound when another prior
     * also depends on the heights.
     */
    public static final class NodeHeightPlan {
        public final String fixedReason;
        public final String treePriorGradient; // element of the tree prior's node-height gradient, or null
        public final String treePriorLikelihood; // element of the tree prior's likelihood
        public final String treePriorId;
        public final boolean numericPrior;
        public final List<Parameter> ctmcPriors;

        private NodeHeightPlan(String fixedReason, String treePriorGradient, String treePriorLikelihood,
                               String treePriorId, boolean numericPrior, List<Parameter> ctmcPriors) {
            this.fixedReason = fixedReason;
            this.treePriorGradient = treePriorGradient;
            this.treePriorLikelihood = treePriorLikelihood;
            this.treePriorId = treePriorId;
            this.numericPrior = numericPrior;
            this.ctmcPriors = ctmcPriors;
        }

        static NodeHeightPlan fixed(String reason) {
            return new NodeHeightPlan(reason, null, null, null, false, new ArrayList<Parameter>());
        }

        public boolean isOptimized() {
            return fixedReason == null;
        }

        public String getDescription() {
            if (fixedReason != null) {
                return "held fixed: " + fixedReason;
            }
            return "analytic (<nodeHeightGradient>" + (numericPrior ? ", numeric prior" :
                    (treePriorGradient != null ? ", <" + treePriorGradient + ">" : "") +
                            (ctmcPriors.isEmpty() ? "" : ", <ctmcScalePriorNodeHeightGradient>")) + ")";
        }
    }

    public static NodeHeightPlan nodeHeightPlan(TreeLimeOptions options, PartitionTreeModel tree) {

        TipDateSamplingComponentOptions tipDates = (TipDateSamplingComponentOptions)
                options.getComponentOptions(TipDateSamplingComponentOptions.class);
        if (tipDates != null && tipDates.tipDateSamplingType != TipDateSamplingType.NO_SAMPLING) {
            return NodeHeightPlan.fixed("gradients with sampled tip dates are not implemented");
        }

        PartitionTreePrior prior = tree.getPartitionTreePrior();
        String treePrefix = tree.getPrefix();
        String gradient = null;
        String likelihood = null;
        String id = null;
        boolean numericPrior = false;

        switch (prior.getNodeHeightPrior()) {
            case CONSTANT:
            case EXPONENTIAL:
            case LOGISTIC:
            case EXPANSION:
                gradient = "coalescentGradient";
                likelihood = "coalescentLikelihood";
                id = treePrefix + "coalescent";
                break;
            case SKYGRID:
            case SKYGRID_HMC:
                gradient = "gmrfSkyrideGradient";
                likelihood = "gmrfSkyGridLikelihood";
                id = prior.getPrefix() + "skygrid";
                // a Skygrid shared by several trees differentiates all their heights together
                numericPrior = treesWithPrior(options, prior) > 1;
                break;
            case GMRF_SKYRIDE:
                gradient = "gmrfSkyrideGradient";
                likelihood = "gmrfSkyrideLikelihood";
                id = treePrefix + "skyride";
                break;
            case YULE:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                gradient = "speciationLikelihoodGradient";
                likelihood = "speciationLikelihood";
                id = treePrefix + "speciation";
                break;
            case YULE_CALIBRATION:
                return NodeHeightPlan.fixed("gradients with calibration priors are not implemented");
            case BIRTH_DEATH_SERIAL_SAMPLING:
                return NodeHeightPlan.fixed("node-height gradients of the serially sampled birth-death model are " +
                        "not implemented");
            default:
                numericPrior = true;
        }
        if (prior.getSubtreeTaxonSet() != null) {
            numericPrior = true;
        }

        List<Parameter> ctmcPriors = new ArrayList<Parameter>();
        for (Parameter parameter : options.selectParameters()) {
            PriorType type = parameter.priorType;
            if ((parameter.isNodeHeight || parameter.isStatistic) && type != PriorType.NONE_TREE_PRIOR &&
                    type != PriorType.NONE_STATISTIC && type != PriorType.UNDEFINED) {
                // a prior on the root height or a tMRCA also depends on the heights
                numericPrior = true;
            }
            if (type == PriorType.CTMC_RATE_REFERENCE_PRIOR) {
                PartitionTreeModel rateTree = ctmcTree(options, parameter);
                if (rateTree == null) {
                    numericPrior = true;
                } else if (rateTree == tree) {
                    ctmcPriors.add(parameter);
                }
            }
        }
        if (numericPrior) {
            return new NodeHeightPlan(null, null, null, null, true, new ArrayList<Parameter>());
        }
        return new NodeHeightPlan(null, gradient, likelihood, id, false, ctmcPriors);
    }

    private static int treesWithPrior(TreeLimeOptions options, PartitionTreePrior prior) {
        int count = 0;
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            if (tree.getPartitionTreePrior() == prior) {
                count++;
            }
        }
        return count;
    }

    /**
     * @return the tree whose length a CTMC rate reference prior uses, as ParameterPriorGenerator finds it, or null
     * when that is not a plain clock rate
     */
    private static PartitionTreeModel ctmcTree(TreeLimeOptions options, Parameter parameter) {
        if (parameter.getTaxonSet() != null) {
            return null;
        }
        PartitionOptions owner = findOwner(options, parameter);
        if (!(owner instanceof PartitionClockModel)) {
            return null;
        }
        PartitionClockModel clock = (PartitionClockModel) owner;
        if (clock.getClockType() == ClockType.FIXED_LOCAL_CLOCK || clock.performModelAveraging()) {
            return null;
        }
        for (AbstractPartitionData partition : options.getDataPartitions(clock)) {
            return partition.getPartitionTreeModel();
        }
        return null;
    }

}
