/*
 * NewFullyConjugateTraitLikelihood.java
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

package dr.app.beagle.evomodel.newtreelikelihood;

import dr.app.beagle.evomodel.parsers.BeagleTreeLikelihoodParser;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ThreadAwareLikelihood;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.Reportable;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
@Deprecated // replaced by TreeDataLikelihood
public class NewFullyConjugateTraitLikelihood extends NewAbstractLikelihoodOnTree implements ThreadAwareLikelihood,
        ConjugateWishartStatisticsProvider, Reportable {

    public NewFullyConjugateTraitLikelihood(String traitName,
                                            TreeModel treeModel,
                                            MultivariateDiffusionModel diffusionModel,
                                            CompoundParameter traitParameter,
                                            Parameter deltaParameter,
                                            List<Integer> missingIndices,
                                            boolean cacheBranches,
                                            boolean scaleByTime,
                                            boolean useTreeLength,
                                            BranchRateModel diffusionRateModel,
                                            List<BranchRateModel> driftRateModels,
                                            BranchRateModel selectionRateModel,
                                            Model samplingDensity,
                                            boolean reportAsMultivariate,
                                            double[] rootPriorMean,
                                            double rootPriorSampleSize,
                                            boolean reciprocalRates,


                                            String name, Map<Set<String>, Parameter> partialsRestrictions) {
        super(BeagleTreeLikelihoodParser.TREE_LIKELIHOOD, // TODO Change
                treeModel, partialsRestrictions);

        this.diffusionModel = diffusionModel;
        this.traitParameter = traitParameter;
        this.diffusionRateModel = diffusionRateModel;
        this.driftRateModels = driftRateModels;
        this.selectionRateModel = selectionRateModel;

        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {
        computeWishartStatistics = true;
        calculateLogLikelihood();
        computeWishartStatistics = false;
        return wishartStatistics;
    }

    final protected MultivariateDiffusionModel diffusionModel;
    final protected CompoundParameter traitParameter;
    final protected BranchRateModel diffusionRateModel;
    final protected List<BranchRateModel> driftRateModels;
    final protected BranchRateModel selectionRateModel;

    protected boolean computeWishartStatistics = false;
    protected WishartSufficientStatistics wishartStatistics;

    // Fully-conjugate-specific
    final double[] rootPriorMean;
    final double rootPriorSampleSize;


}
