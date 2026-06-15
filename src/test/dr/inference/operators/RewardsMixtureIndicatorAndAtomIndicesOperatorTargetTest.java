/*
 * RewardsMixtureIndicatorAndAtomIndicesOperatorTargetTest.java
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.operators;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.BeagleRewardDependentCtmcEdgeEvidenceProvider;
import dr.inference.operators.RewardDependentCtmcEdgeEvidenceProvider;
import dr.inference.operators.RewardsMixtureBranchResamplingHelper;
import dr.inference.operators.RewardsMixtureClusterResamplingHelper;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.math.MathUtils;
import test.dr.math.MathTestCase;

import java.lang.reflect.Method;
import java.util.Arrays;

/**
 * Exact target regressions for the reward-mixture operator when dependent CTMC
 * likelihoods contribute branch-local evidence.
 *
 * @author Filippo Monti
 */
public class RewardsMixtureIndicatorAndAtomIndicesOperatorTargetTest extends MathTestCase {

    private static final double TOL = 1.0e-7;

    public void testSingleBranchWeightsMatchFullTargetWithOneDependentCtmc() throws Exception {
        final Fixture fixture = createFixture(1);
        assertSingleBranchWeightsMatchFullTarget(fixture, firstCherryExternalNode(fixture.tree));
    }

    public void testSingleBranchWeightsSumTwoDependentCtmcLikelihoods() throws Exception {
        final Fixture fixture = createFixture(2);
        assertSingleBranchWeightsMatchFullTarget(fixture, firstCherryExternalNode(fixture.tree));
    }

    public void testClusterProposalTargetsMatchFullLikelihoodWithDependentCtmc() {
        MathUtils.setSeed(20260611);

        final Fixture fixture = createFixture(1);
        final RewardDependentCtmcEdgeEvidenceProvider[] providers = createDependentEvidenceProviders(fixture);
        final RewardsMixtureClusterResamplingHelper helper = new RewardsMixtureClusterResamplingHelper(
                fixture.indicator,
                fixture.atomIndices,
                fixture.rewardsAwareBranchModel,
                fixture.independentLikelihood,
                fixture.independentDelegate,
                providers,
                2,
                1.0,
                new RewardsMixtureClusterResamplingHelper.ExactLogTargetEvaluator() {
                    @Override
                    public double computeCurrentLogTarget(final int[] clusterBranchNodeNumbers,
                                                          final int clusterCount) {
                        return currentFullLogTarget(fixture);
                    }
                },
                new RewardsMixtureClusterResamplingHelper.LikelihoodStateRefresher() {
                    @Override
                    public void refreshCurrentState() {
                        refreshLikelihoodMessages(fixture, providers);
                    }
                }
        );

        final int seedBranchNodeNumber = firstCherryExternalNode(fixture.tree).getNumber();
        refreshLikelihoodMessages(fixture, providers);
        final double exactOld = currentFullLogTarget(fixture);

        final RewardsMixtureClusterResamplingHelper.ClusterProposal proposal =
                helper.proposeClusterMove(seedBranchNodeNumber);

        assertTrue("Expected at least one branch in the proposed cluster", proposal.clusterCount >= 1);
        assertTrue("Cluster exceeds requested maximum size", proposal.clusterCount <= 2);
        assertFinite(proposal.logTargetOld);
        assertFinite(proposal.logTargetNew);
        assertFinite(proposal.logQForward);
        assertFinite(proposal.logQReverse);
        assertEquals(exactOld, proposal.logTargetOld, TOL);

        final double exactNew = recomputeFullLogTarget(fixture);
        assertEquals(exactNew, proposal.logTargetNew, TOL);
        assertProposalStateIsApplied(fixture, proposal);

        proposal.restoreOldState();
        proposal.fireParameterEvents();
        markLikelihoodsDirty(fixture);
    }

    private static void assertSingleBranchWeightsMatchFullTarget(final Fixture fixture,
                                                                 final NodeRef branch) throws Exception {
        final int branchNodeNumber = branch.getNumber();
        final int branchIndex = fixture.rewardsAwareBranchModel.getParameterIndexForNode(branchNodeNumber);
        final double[] oldIndicators = copyParameterValues(fixture.indicator);
        final double[] oldAtoms = copyParameterValues(fixture.atomIndices);

        try {
            setBranchContinuous(fixture, branchIndex);
            markLikelihoodsDirty(fixture);

            final RewardsMixtureIndicatorAndAtomIndicesOperator operator =
                    new RewardsMixtureIndicatorAndAtomIndicesOperator(
                            fixture.indicator,
                            fixture.atomIndices,
                            fixture.rewardsAwareBranchModel,
                            fixture.independentLikelihood,
                            fixture.dependentLikelihoods,
                            1.0,
                            false,
                            false,
                            1,
                            0.5,
                            1.0
                    );

            final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                    computeBranchWeights(fixture, operator, branchNodeNumber);

            final double logTargetCts = recomputeFullLogTarget(fixture);
            final int stateCount = fixture.rewardsAwareBranchModel.getStateCount();
            final double[] logAtomicTargets = new double[stateCount];

            for (int state = 0; state < stateCount; state++) {
                setBranchAtomic(fixture, branchIndex, state);
                logAtomicTargets[state] = recomputeFullLogTargetAllowingNegativeInfinity(fixture);

                assertLogEquals(
                        "Full target ratio must match local weight ratio for atom state " + state,
                        logAtomicTargets[state] - logTargetCts,
                        weights.logAtomicWeights[state] - weights.logCtsWeight
                );
            }

            assertLogEquals(
                    "Summed atomic full target must match summed atomic local weight",
                    logSum(logAtomicTargets) - logTargetCts,
                    weights.logAtomicTotalWeight - weights.logCtsWeight
            );
        } finally {
            restoreParameterValues(fixture.indicator, oldIndicators);
            restoreParameterValues(fixture.atomIndices, oldAtoms);
            markLikelihoodsDirty(fixture);
        }
    }

    private static RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeights(
            final Fixture fixture,
            final RewardsMixtureIndicatorAndAtomIndicesOperator operator,
            final int branchNodeNumber
    ) throws Exception {
        markLikelihoodsDirty(fixture);

        final Method refresh = RewardsMixtureIndicatorAndAtomIndicesOperator.class
                .getDeclaredMethod("refreshLikelihoodMessages");
        refresh.setAccessible(true);
        refresh.invoke(operator);

        final Method compute = RewardsMixtureIndicatorAndAtomIndicesOperator.class
                .getDeclaredMethod("computeBranchWeights", int.class);
        compute.setAccessible(true);
        return (RewardsMixtureBranchResamplingHelper.BranchWeights) compute.invoke(operator, branchNodeNumber);
    }

    private static Fixture createFixture(final int dependentCount) {
        final TreeModel tree = createThreeTipTree();
        final SubstitutionModel independentSubstitutionModel =
                createNucleotideSubstitutionModel("rewardMixtureTargetIndependentKappa", 2.0);
        final SitePatterns independentPatterns = createSitePatterns("A", "C", "G");
        final GammaSiteRateModel independentSiteRateModel =
                new GammaSiteRateModel("rewardMixtureTargetIndependentRateModel");

        final Parameter ctsRewards =
                new Parameter.Default("targetRewardCts", new double[]{0.50, 0.50, 0.50, 0.50});
        final Parameter indicator =
                new Parameter.Default("targetRewardIndicator", new double[]{0.0, 0.0, 0.0, 0.0});
        final Parameter atomIndices =
                new Parameter.Default("targetRewardAtomIndices", new double[]{0.0, 1.0, 2.0, 3.0});
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("targetRewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("targetRewardRatesInternal", new double[0]),
                new Parameter.Default("targetRewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
        final RewardsAwareMixtureBranchRates rewardBranchRates = new RewardsAwareMixtureBranchRates(
                tree,
                ctsRewards,
                indicator,
                atomIndices,
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                independentSubstitutionModel,
                rewardRates,
                indicator,
                rewardBranchRates,
                atomIndices,
                false
        );

        final DiscreteDataLikelihoodDelegate independentDelegate = new DiscreteDataLikelihoodDelegate(
                tree,
                independentPatterns,
                rewardsAwareBranchModel,
                independentSiteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, false, false, false),
                new RewardsAwarePartialsRepresentation(rewardsAwareBranchModel),
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                null,
                null
        );
        final TreeDataLikelihood independentLikelihood =
                new TreeDataLikelihood(independentDelegate, tree, new DefaultBranchRateModel());

        final TreeDataLikelihood[] dependentLikelihoods = new TreeDataLikelihood[dependentCount];
        for (int i = 0; i < dependentCount; i++) {
            dependentLikelihoods[i] = createDependentLikelihood(tree, rewardBranchRates, i);
        }

        return new Fixture(
                tree,
                rewardsAwareBranchModel,
                independentDelegate,
                independentLikelihood,
                dependentLikelihoods,
                indicator,
                atomIndices
        );
    }

    private static TreeDataLikelihood createDependentLikelihood(final TreeModel tree,
                                                                final RewardsAwareMixtureBranchRates rewardBranchRates,
                                                                final int index) {
        final String[] taxonStates = index % 2 == 0
                ? new String[]{"A", "C", "G"}
                : new String[]{"T", "A", "C"};
        final double kappa = index % 2 == 0 ? 3.0 : 1.7;
        final SubstitutionModel substitutionModel =
                createNucleotideSubstitutionModel("rewardMixtureDependentKappa" + index, kappa);
        final GammaSiteRateModel siteRateModel =
                new GammaSiteRateModel("rewardMixtureDependentRateModel" + index);
        final BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(
                tree,
                createSitePatterns(taxonStates[0], taxonStates[1], taxonStates[2]),
                new HomogeneousBranchModel(substitutionModel),
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, true, false, false, false)
        );
        return new TreeDataLikelihood(delegate, tree, rewardBranchRates);
    }

    private static SitePatterns createSitePatterns(final String a, final String b, final String c) {
        final DataType dataType = Nucleotides.INSTANCE;
        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        alignment.addSequence(sequence("a", a, dataType));
        alignment.addSequence(sequence("b", b, dataType));
        alignment.addSequence(sequence("c", c, dataType));
        return new SitePatterns(alignment, null, 0, -1, 1, true);
    }

    private static Sequence sequence(final String taxonName, final String sequenceString, final DataType dataType) {
        final Sequence sequence = new Sequence(new Taxon(taxonName), sequenceString);
        sequence.setDataType(dataType);
        return sequence;
    }

    private static TreeModel createThreeTipTree() {
        final SimpleNode left = new SimpleNode();
        left.setTaxon(new Taxon("a"));
        left.setHeight(0.0);

        final SimpleNode right = new SimpleNode();
        right.setTaxon(new Taxon("b"));
        right.setHeight(0.0);

        final SimpleNode internal = new SimpleNode();
        internal.setHeight(0.4);
        internal.addChild(left);
        internal.addChild(right);

        final SimpleNode outgroup = new SimpleNode();
        outgroup.setTaxon(new Taxon("c"));
        outgroup.setHeight(0.0);

        final SimpleNode root = new SimpleNode();
        root.setHeight(1.0);
        root.addChild(internal);
        root.addChild(outgroup);

        final Tree tree = new SimpleTree(root);
        return new DefaultTreeModel("rewardMixtureTargetTree", tree);
    }

    private static SubstitutionModel createNucleotideSubstitutionModel(final String id, final double kappa) {
        final DataType dataType = Nucleotides.INSTANCE;
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.25, 0.25, 0.25, 0.25});
        return new HKY(new Parameter.Default(id, kappa), frequencyModel);
    }

    private static NodeRef firstCherryExternalNode(final Tree tree) {
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            final NodeRef node = tree.getExternalNode(i);
            final NodeRef parent = tree.getParent(node);
            if (parent != null && !tree.isRoot(parent)) {
                return node;
            }
        }
        throw new IllegalArgumentException("Tree has no external branch below a non-root internal node");
    }

    private static RewardDependentCtmcEdgeEvidenceProvider[] createDependentEvidenceProviders(final Fixture fixture) {
        final RewardDependentCtmcEdgeEvidenceProvider[] providers =
                new RewardDependentCtmcEdgeEvidenceProvider[fixture.dependentLikelihoods.length];
        for (int i = 0; i < providers.length; i++) {
            providers[i] = new BeagleRewardDependentCtmcEdgeEvidenceProvider(fixture.dependentLikelihoods[i]);
        }
        return providers;
    }

    private static void refreshLikelihoodMessages(final Fixture fixture,
                                                  final RewardDependentCtmcEdgeEvidenceProvider[] providers) {
        markLikelihoodsDirty(fixture);
        fixture.independentDelegate.updatePostOrdersFromTreeDataLikelihood(fixture.independentLikelihood);
        fixture.independentDelegate.ensurePreOrderComputed();
        for (RewardDependentCtmcEdgeEvidenceProvider provider : providers) {
            provider.prepare();
        }
    }

    private static double recomputeFullLogTarget(final Fixture fixture) {
        markLikelihoodsDirty(fixture);
        return currentFullLogTarget(fixture);
    }

    private static double recomputeFullLogTargetAllowingNegativeInfinity(final Fixture fixture) {
        markLikelihoodsDirty(fixture);
        return currentFullLogTarget(fixture, false);
    }

    private static double currentFullLogTarget(final Fixture fixture) {
        return currentFullLogTarget(fixture, true);
    }

    private static double currentFullLogTarget(final Fixture fixture, final boolean requireFinite) {
        double logTarget = fixture.independentLikelihood.getLogLikelihood();
        for (TreeDataLikelihood dependentLikelihood : fixture.dependentLikelihoods) {
            logTarget += dependentLikelihood.getLogLikelihood();
        }
        if (requireFinite) {
            assertFinite(logTarget);
        } else {
            assertTrue("Expected non-NaN target but found " + logTarget, !Double.isNaN(logTarget));
            assertTrue("Unexpected positive infinite target", logTarget != Double.POSITIVE_INFINITY);
        }
        return logTarget;
    }

    private static void markLikelihoodsDirty(final Fixture fixture) {
        fixture.independentLikelihood.makeDirty();
        for (TreeDataLikelihood dependentLikelihood : fixture.dependentLikelihoods) {
            dependentLikelihood.makeDirty();
        }
    }

    private static void setBranchContinuous(final Fixture fixture, final int branchIndex) {
        fixture.indicator.setParameterValue(branchIndex, 0.0);
    }

    private static void setBranchAtomic(final Fixture fixture, final int branchIndex, final int state) {
        fixture.atomIndices.setParameterValue(branchIndex, state);
        fixture.indicator.setParameterValue(branchIndex, 1.0);
    }

    private static void assertProposalStateIsApplied(final Fixture fixture,
                                                     final RewardsMixtureClusterResamplingHelper.ClusterProposal proposal) {
        for (int i = 0; i < proposal.clusterCount; i++) {
            final int branchIndex = proposal.parameterIndices[i];
            assertEquals(proposal.newIndicators[i],
                    (int) Math.round(fixture.indicator.getParameterValue(branchIndex)));
            if (proposal.newIndicators[i] == 1) {
                assertEquals(proposal.newAtoms[i],
                        (int) Math.round(fixture.atomIndices.getParameterValue(branchIndex)));
            }
        }
    }

    private static double[] copyParameterValues(final Parameter parameter) {
        final double[] values = new double[parameter.getDimension()];
        for (int i = 0; i < values.length; i++) {
            values[i] = parameter.getParameterValue(i);
        }
        return values;
    }

    private static void restoreParameterValues(final Parameter parameter, final double[] values) {
        for (int i = 0; i < values.length; i++) {
            parameter.setParameterValue(i, values[i]);
        }
    }

    private static double logSum(final double[] x) {
        double acc = Double.NEGATIVE_INFINITY;
        for (double v : x) {
            acc = RewardsMixtureBranchResamplingHelper.logAdd(acc, v);
        }
        return acc;
    }

    private static void assertLogEquals(final String message,
                                        final double expected,
                                        final double actual) {
        if (Double.isInfinite(expected) || Double.isInfinite(actual)) {
            assertEquals(message, expected, actual);
        } else {
            assertEquals(message, expected, actual, TOL);
        }
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static final class Fixture {
        final TreeModel tree;
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final DiscreteDataLikelihoodDelegate independentDelegate;
        final TreeDataLikelihood independentLikelihood;
        final TreeDataLikelihood[] dependentLikelihoods;
        final Parameter indicator;
        final Parameter atomIndices;

        private Fixture(final TreeModel tree,
                        final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final DiscreteDataLikelihoodDelegate independentDelegate,
                        final TreeDataLikelihood independentLikelihood,
                        final TreeDataLikelihood[] dependentLikelihoods,
                        final Parameter indicator,
                        final Parameter atomIndices) {
            this.tree = tree;
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.independentDelegate = independentDelegate;
            this.independentLikelihood = independentLikelihood;
            this.dependentLikelihoods = Arrays.copyOf(dependentLikelihoods, dependentLikelihoods.length);
            this.indicator = indicator;
            this.atomIndices = atomIndices;
        }
    }
}
