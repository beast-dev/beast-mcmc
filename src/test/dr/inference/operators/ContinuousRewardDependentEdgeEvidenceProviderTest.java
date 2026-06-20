/*
 * ContinuousRewardDependentEdgeEvidenceProviderTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
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
 */

package test.dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.HomogeneousDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Parameter;
import dr.inference.operators.BranchLocalContinuousRewardDependentEdgeEvidenceProvider;
import dr.inference.operators.ContinuousRewardDependentEdgeEvidenceProvider;
import test.dr.evomodel.treedatalikelihood.continuous.ContinuousTraitTest;

import java.util.ArrayList;
import java.util.List;

/**
 * Exact continuous edge evidence should match a direct full likelihood
 * evaluation under the same candidate raw reward and restore the current state.
 */
public class ContinuousRewardDependentEdgeEvidenceProviderTest extends ContinuousTraitTest {

    private static final double TOL = 1.0e-10;
    private static final double DELTA_TOL = 1.0e-8;

    public ContinuousRewardDependentEdgeEvidenceProviderTest(String name) {
        super(name);
    }

    public void testProviderMatchesManualContinuousLikelihoodEvaluation() {
        final RewardsAwareMixtureBranchRates rewardBranchRates = newRewardBranchRates();
        final Parameter ctsRewards = rewardBranchRates.getRateParameter();
        final Parameter indicator = rewardBranchRates.getIndicator();

        final DiffusionProcessDelegate diffusionProcessDelegate =
                new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                rateTransformation,
                rewardBranchRates,
                true);
        final TreeDataLikelihood dataLikelihood =
                new TreeDataLikelihood(likelihoodDelegate, treeModel, rewardBranchRates);
        final ContinuousRewardDependentEdgeEvidenceProvider provider =
                new ContinuousRewardDependentEdgeEvidenceProvider(dataLikelihood);

        final NodeRef branch = firstNonRootNode();
        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(branch);
        final double oldReward = ctsRewards.getParameterValue(parameterIndex);
        final double oldIndicator = indicator.getParameterValue(parameterIndex);
        final double candidateReward = 0.73;

        final double expected = manualLogLikelihood(
                dataLikelihood, ctsRewards, indicator, parameterIndex, candidateReward);
        final double actual = provider.logEvidence(branch.getNumber(), candidateReward);

        assertEquals(expected, actual, TOL);
        assertEquals(oldReward, ctsRewards.getParameterValue(parameterIndex), 0.0);
        assertEquals(oldIndicator, indicator.getParameterValue(parameterIndex), 0.0);
    }

    public void testBranchLocalProviderMatchesExactBrownianDeltas() {
        final RewardsAwareMixtureBranchRates rewardBranchRates = newRewardBranchRates();
        final TreeDataLikelihood localLikelihood = newBrownianLikelihood(rewardBranchRates);
        final TreeDataLikelihood exactLikelihood = newBrownianLikelihood(rewardBranchRates);

        assertLocalProviderMatchesExactDeltas(localLikelihood, exactLikelihood);
    }

    public void testBranchLocalProviderMatchesExactDiagonalOUDeltas() {
        final RewardsAwareMixtureBranchRates rewardBranchRates = newRewardBranchRates();
        final TreeDataLikelihood localLikelihood = newOULikelihood(rewardBranchRates);
        final TreeDataLikelihood exactLikelihood = newOULikelihood(rewardBranchRates);

        assertLocalProviderMatchesExactDeltas(localLikelihood, exactLikelihood);
    }

    private NodeRef firstNonRootNode() {
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            final NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
                return node;
            }
        }
        throw new IllegalStateException("Tree has no non-root branch");
    }

    private void assertLocalProviderMatchesExactDeltas(final TreeDataLikelihood localLikelihood,
                                                       final TreeDataLikelihood exactLikelihood) {
        final BranchLocalContinuousRewardDependentEdgeEvidenceProvider localProvider =
                new BranchLocalContinuousRewardDependentEdgeEvidenceProvider(localLikelihood);
        final ContinuousRewardDependentEdgeEvidenceProvider exactProvider =
                new ContinuousRewardDependentEdgeEvidenceProvider(exactLikelihood);

        final double referenceReward = 0.55;
        final double[] candidateRewards = new double[]{0.25, 0.73, 1.15};

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            final NodeRef node = treeModel.getNode(i);
            if (treeModel.isRoot(node)) {
                continue;
            }

            localProvider.prepare();
            final double localReference = localProvider.logEvidence(node.getNumber(), referenceReward);

            for (double candidateReward : candidateRewards) {
                final double localDelta =
                        localProvider.logEvidence(node.getNumber(), candidateReward) - localReference;
                final double exactDelta =
                        exactProvider.logEvidence(node.getNumber(), candidateReward) -
                                exactProvider.logEvidence(node.getNumber(), referenceReward);
                assertEquals("node " + node.getNumber() + " candidate " + candidateReward,
                        exactDelta, localDelta, DELTA_TOL);
            }
        }
    }

    private TreeDataLikelihood newBrownianLikelihood(final RewardsAwareMixtureBranchRates rewardBranchRates) {
        final DiffusionProcessDelegate diffusionProcessDelegate =
                new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                rateTransformation,
                rewardBranchRates,
                true);
        return new TreeDataLikelihood(likelihoodDelegate, treeModel, rewardBranchRates);
    }

    private TreeDataLikelihood newOULikelihood(final RewardsAwareMixtureBranchRates rewardBranchRates) {
        final List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("ou.optimum.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("ou.optimum.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("ou.optimum.3", new double[]{-2.0})));

        final DiagonalMatrix strengthOfSelectionMatrixParam =
                new DiagonalMatrix(new Parameter.Default(new double[]{0.2, 0.5, 0.8}));

        final DiffusionProcessDelegate diffusionProcessDelegate =
                new OUDiffusionModelDelegate(treeModel, diffusionModel,
                        optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                rateTransformation,
                rewardBranchRates,
                true);
        return new TreeDataLikelihood(likelihoodDelegate, treeModel, rewardBranchRates);
    }

    private RewardsAwareMixtureBranchRates newRewardBranchRates() {
        final int branchCount = treeModel.getNodeCount() - 1;
        final Parameter ctsRewards = new Parameter.Default("ctsRewards", fill(branchCount, 0.50));
        final Parameter indicator = new Parameter.Default("indicator", fill(branchCount, 1.0));
        final Parameter atomIndices = new Parameter.Default("atomIndices", alternatingAtoms(branchCount));
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("rewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("rewardRatesInternal", new double[0]),
                new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );

        return new RewardsAwareMixtureBranchRates(
                treeModel,
                ctsRewards,
                indicator,
                atomIndices,
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );
    }

    private static double manualLogLikelihood(final TreeDataLikelihood likelihood,
                                              final Parameter rewards,
                                              final Parameter indicator,
                                              final int parameterIndex,
                                              final double candidateReward) {
        final double oldReward = rewards.getParameterValue(parameterIndex);
        final double oldIndicator = indicator.getParameterValue(parameterIndex);
        try {
            rewards.setParameterValueQuietly(parameterIndex, candidateReward);
            indicator.setParameterValueQuietly(parameterIndex, 0.0);
            likelihood.makeDirty();
            return likelihood.getLogLikelihood();
        } finally {
            rewards.setParameterValueQuietly(parameterIndex, oldReward);
            indicator.setParameterValueQuietly(parameterIndex, oldIndicator);
            likelihood.makeDirty();
        }
    }

    private static double[] fill(final int n, final double value) {
        final double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = value;
        }
        return x;
    }

    private static double[] alternatingAtoms(final int n) {
        final double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i % 4;
        }
        return x;
    }
}
