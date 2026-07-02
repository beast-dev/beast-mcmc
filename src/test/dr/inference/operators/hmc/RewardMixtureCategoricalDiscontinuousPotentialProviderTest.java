/*
 * RewardMixtureCategoricalDiscontinuousPotentialProviderTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

package test.dr.inference.operators.hmc;

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
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRates;
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
import dr.inference.hmc.RewardMixtureCategoricalDiscontinuousPotentialProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.BeagleRewardDependentCtmcEdgeEvidenceProvider;
import dr.inference.operators.RewardDependentEdgeEvidenceProvider;
import dr.inference.operators.RewardsMixtureBranchWeightProvider;
import dr.inference.operators.hmc.DiscontinuousHamiltonianMonteCarloOperator;
import dr.math.MathUtils;
import test.dr.math.MathTestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Target-ratio regressions for categorical reward-mixture DHMC potential jumps.
 */
public class RewardMixtureCategoricalDiscontinuousPotentialProviderTest extends MathTestCase {

    private static final double TOL = 1.0e-7;

    public void testCategoricalBranchRateDecoding() {
        final Fixture fixture = createFixture();
        final int branchIndex = 0;
        final NodeRef node = fixture.tree.getNode(
                fixture.branchRates.getNodeNumberFromParameterIndex(branchIndex));

        setCategory(fixture, branchIndex, 0);
        assertEquals(fixture.ctsRewards.getParameterValue(branchIndex),
                fixture.branchRates.getUntransformedBranchRate(fixture.tree, node), TOL);
        assertFalse(fixture.rewardsAwareBranchModel.isAtomicBranch(node.getNumber()));

        setCategory(fixture, branchIndex, 3);
        assertEquals(0.60, fixture.branchRates.getUntransformedBranchRate(fixture.tree, node), TOL);
        assertEquals(0.0, fixture.branchRates.getBranchRateDifferential(fixture.tree, node), TOL);
        assertTrue(fixture.rewardsAwareBranchModel.isAtomicBranch(node.getNumber()));
        assertEquals(2, fixture.rewardsAwareBranchModel.getAtomicBranchState(node.getNumber()));
    }

    public void testPotentialDifferenceMatchesFullTargetRatiosFromContinuousState() {
        final Fixture fixture = createFixture();
        final NodeRef branch = firstCherryExternalNode(fixture.tree);
        final int branchIndex = fixture.rewardsAwareBranchModel.getParameterIndexForNode(branch.getNumber());
        final RewardMixtureCategoricalDiscontinuousPotentialProvider provider = createProvider(fixture);

        final double currentValue = categoryValue(0);

        for (int category = 0; category < provider.getCategoryCount(); category++) {
            setCategory(fixture, branchIndex, 0);
            markLikelihoodDirty(fixture);
            final double currentTarget = currentFullLogTarget(fixture);

            final double proposedValue = categoryValue(category);
            final double deltaU =
                    provider.getPotentialDifference(branchIndex, currentValue, proposedValue);

            setCategory(fixture, branchIndex, category);
            markLikelihoodDirty(fixture);
            final double proposedTarget = currentFullLogTargetAllowingNegativeInfinity(fixture);

            assertLogEquals(
                    "Potential jump must match full target ratio for category " + category,
                    currentTarget - proposedTarget,
                    deltaU
            );
        }
    }

    public void testContinuousCandidateUsesLatentContinuousRewardFromAtomicState() {
        final Fixture fixture = createFixture();
        final NodeRef branch = firstCherryExternalNode(fixture.tree);
        final int branchIndex = fixture.rewardsAwareBranchModel.getParameterIndexForNode(branch.getNumber());
        final RewardMixtureCategoricalDiscontinuousPotentialProvider provider = createProvider(fixture);

        fixture.ctsRewards.setParameterValue(branchIndex, 0.55);
        setCategory(fixture, branchIndex, 1);
        markLikelihoodDirty(fixture);

        final double atomicValue = categoryValue(1);
        final double atomicTarget = currentFullLogTarget(fixture);
        final double deltaU = provider.getPotentialDifference(branchIndex, atomicValue, categoryValue(0));

        setCategory(fixture, branchIndex, 0);
        markLikelihoodDirty(fixture);
        final double continuousTarget = currentFullLogTarget(fixture);

        assertLogEquals(
                "Continuous candidate from an atomic state must use the latent continuous reward",
                atomicTarget - continuousTarget,
                deltaU
        );
    }

    public void testLogDensityAfterSingleCoordinateMoveIsConsistentWithPotentialDifference() {
        final Fixture fixture = createFixture();
        final NodeRef branch = firstCherryExternalNode(fixture.tree);
        final int branchIndex = fixture.rewardsAwareBranchModel.getParameterIndexForNode(branch.getNumber());
        final RewardMixtureCategoricalDiscontinuousPotentialProvider provider = createProvider(fixture);

        setCategory(fixture, branchIndex, 0);
        final double currentLogDensity = provider.getLogDensity();
        final double proposedLogDensity =
                provider.getLogDensityAfterSingleCoordinateMove(branchIndex, categoryValue(1));
        final double deltaU =
                provider.getPotentialDifference(branchIndex, categoryValue(0), categoryValue(1));

        assertLogEquals(currentLogDensity - proposedLogDensity, deltaU);
    }

    public void testDiscontinuousHmcOperatorRunsOnRewardAwareCategories() {
        MathUtils.setSeed(20260630);

        final Fixture fixture = createFixture();
        final RewardMixtureCategoricalDiscontinuousPotentialProvider provider = createProvider(fixture);
        final double[] momentumScales = new double[fixture.categories.getDimension()];
        Arrays.fill(momentumScales, 1.0);

        final DiscontinuousHamiltonianMonteCarloOperator operator =
                new DiscontinuousHamiltonianMonteCarloOperator(
                        provider,
                        momentumScales,
                        1.0,
                        3,
                        1.0);

        assertFinite(provider.getLogDensity());

        for (int i = 0; i < 20; i++) {
            assertEquals(0.0, operator.doOperation(), TOL);
            markLikelihoodDirty(fixture);
            assertFinite(currentFullLogTarget(fixture));

            for (int branch = 0; branch < fixture.categories.getDimension(); branch++) {
                final int category = provider.getRewardMixtureCategory(branch);
                assertTrue("DHMC category out of range: " + category,
                        category >= 0 && category < provider.getCategoryCount());
            }
        }

        assertTrue("DHMC should have crossed or reflected at least once",
                operator.getCrossingCount() + operator.getReflectionCount() > 0);
    }

    public void testDiscontinuousHmcWithDependentCtmcDiagnosticsMatchesExactDeltas() throws Exception {
        MathUtils.setSeed(20260630);

        final Fixture fixture = createFixtureWithDependentCtmc();
        final File diagnosticFile =
                File.createTempFile("categorical_dhmc_dependent_ctmc_edge_evidence", ".tsv");
        diagnosticFile.deleteOnExit();

        final BeagleRewardDependentCtmcEdgeEvidenceProvider.Diagnostics diagnostics =
                BeagleRewardDependentCtmcEdgeEvidenceProvider.Diagnostics.create(
                        true,
                        diagnosticFile.getAbsolutePath(),
                        true,
                        true,
                        Integer.MAX_VALUE,
                        128);
        final RewardMixtureCategoricalDiscontinuousPotentialProvider provider =
                createProvider(fixture, diagnostics);

        final double[] momentumScales = new double[fixture.categories.getDimension()];
        Arrays.fill(momentumScales, 1.0);
        final DiscontinuousHamiltonianMonteCarloOperator operator =
                new DiscontinuousHamiltonianMonteCarloOperator(
                        provider,
                        momentumScales,
                        0.25,
                        2,
                        1.0);

        assertFinite(currentJointLogTarget(fixture));
        assertFinite(provider.getLogDensity());

        for (int i = 0; i < 4; i++) {
            assertEquals(0.0, operator.doOperation(), TOL);
            markAllLikelihoodsDirty(fixture);
            assertFinite(currentJointLogTarget(fixture));
            assertFinite(provider.getLogDensity());

            for (int branch = 0; branch < fixture.categories.getDimension(); branch++) {
                final int category = provider.getRewardMixtureCategory(branch);
                assertTrue("DHMC category out of range: " + category,
                        category >= 0 && category < provider.getCategoryCount());
            }
        }

        assertDiagnosticDeltasMatch(diagnosticFile);
    }

    private static RewardMixtureCategoricalDiscontinuousPotentialProvider createProvider(final Fixture fixture) {
        return new RewardMixtureCategoricalDiscontinuousPotentialProvider(
                fixture.categories,
                fixture.categoryCuts,
                fixture.rewardsAwareBranchModel,
                fixture.independentLikelihood,
                null,
                null
        );
    }

    private static RewardMixtureCategoricalDiscontinuousPotentialProvider createProvider(
            final Fixture fixture,
            final BeagleRewardDependentCtmcEdgeEvidenceProvider.Diagnostics dependentCtmcDiagnostics) {
        if (fixture.dependentLikelihood == null) {
            throw new IllegalArgumentException("Fixture has no dependent CTMC likelihood");
        }
        final RewardDependentEdgeEvidenceProvider[] dependentEvidenceProviders =
                new RewardDependentEdgeEvidenceProvider[] {
                        new BeagleRewardDependentCtmcEdgeEvidenceProvider(
                                fixture.dependentLikelihood,
                                dependentCtmcDiagnostics)
                };
        final RewardsMixtureBranchWeightProvider branchWeightProvider =
                new RewardsMixtureBranchWeightProvider(
                        fixture.rewardsAwareBranchModel,
                        fixture.independentLikelihood,
                        (DiscreteDataLikelihoodDelegate) fixture.independentLikelihood.getDataLikelihoodDelegate(),
                        dependentEvidenceProviders);
        return new RewardMixtureCategoricalDiscontinuousPotentialProvider(
                fixture.categories,
                fixture.categoryCuts,
                branchWeightProvider);
    }

    private static Fixture createFixtureWithDependentCtmc() {
        final Fixture fixture = createFixture();
        return fixture.withDependentLikelihood(
                createDependentCtmcLikelihood(fixture.tree, fixture.branchRates));
    }

    private static Fixture createFixture() {
        final TreeModel tree = createThreeTipTree();
        final SubstitutionModel independentSubstitutionModel =
                createNucleotideSubstitutionModel("categoricalRewardMixtureKappa", 2.0);
        final SitePatterns independentPatterns = createSitePatterns("A", "C", "G");
        final GammaSiteRateModel independentSiteRateModel =
                new GammaSiteRateModel("categoricalRewardMixtureRateModel");

        final Parameter ctsRewards =
                new Parameter.Default("categoricalRewardCts", new double[]{0.50, 0.45, 0.55, 0.65});
        final Parameter categories =
                new Parameter.Default("categoricalRewardState", new double[]{0.5, 0.5, 0.5, 0.5});
        final Parameter categoryCuts =
                new Parameter.Default("categoricalRewardCuts", new double[]{0.0, 1.0, 2.0, 3.0, 4.0, 5.0});
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("categoricalRewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("categoricalRewardRatesInternal", new double[0]),
                new Parameter.Default("categoricalRewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );

        final RewardsAwareCategoricalMixtureBranchRates branchRates =
                new RewardsAwareCategoricalMixtureBranchRates(
                        tree,
                        ctsRewards,
                        categories,
                        categoryCuts,
                        rewardRates,
                        new ArbitraryBranchRates.BranchRateTransform.None(),
                        false,
                        TreeParameterModel.Type.WITHOUT_ROOT
                );
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                independentSubstitutionModel,
                rewardRates,
                categories,
                categoryCuts,
                branchRates,
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

        return new Fixture(
                tree,
                ctsRewards,
                categories,
                categoryCuts,
                branchRates,
                rewardsAwareBranchModel,
                independentLikelihood
        );
    }

    private static TreeDataLikelihood createDependentCtmcLikelihood(
            final TreeModel tree,
            final RewardsAwareCategoricalMixtureBranchRates branchRates) {
        final SubstitutionModel substitutionModel =
                createNucleotideSubstitutionModel("categoricalRewardDependentKappa", 2.0);
        final SitePatterns patterns = createSitePatterns("ACGT", "CGTA", "GTAC");
        final GammaSiteRateModel siteRateModel =
                new GammaSiteRateModel("categoricalRewardDependentRateModel");
        final BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(
                tree,
                patterns,
                new HomogeneousBranchModel(substitutionModel),
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, true, false, false, false)
        );
        return new TreeDataLikelihood(delegate, tree, branchRates);
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
        return new DefaultTreeModel("categoricalRewardMixtureTree", tree);
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

    private static void setCategory(final Fixture fixture, final int branchIndex, final int category) {
        fixture.categories.setParameterValue(branchIndex, categoryValue(category));
    }

    private static double categoryValue(final int category) {
        return category + 0.5;
    }

    private static double currentFullLogTarget(final Fixture fixture) {
        final double logTarget = fixture.independentLikelihood.getLogLikelihood();
        assertFinite(logTarget);
        return logTarget;
    }

    private static double currentJointLogTarget(final Fixture fixture) {
        if (fixture.dependentLikelihood == null) {
            return currentFullLogTarget(fixture);
        }
        final double independent = fixture.independentLikelihood.getLogLikelihood();
        final double dependent = fixture.dependentLikelihood.getLogLikelihood();
        assertFinite(independent);
        assertFinite(dependent);
        return independent + dependent;
    }

    private static double currentFullLogTargetAllowingNegativeInfinity(final Fixture fixture) {
        final double logTarget = fixture.independentLikelihood.getLogLikelihood();
        assertTrue("Expected non-NaN target but found " + logTarget, !Double.isNaN(logTarget));
        assertTrue("Unexpected positive infinite target", logTarget != Double.POSITIVE_INFINITY);
        return logTarget;
    }

    private static void markLikelihoodDirty(final Fixture fixture) {
        fixture.independentLikelihood.makeDirty();
    }

    private static void markAllLikelihoodsDirty(final Fixture fixture) {
        markLikelihoodDirty(fixture);
        if (fixture.dependentLikelihood != null) {
            fixture.dependentLikelihood.makeDirty();
        }
    }

    private static void assertDiagnosticDeltasMatch(final File diagnosticFile) throws Exception {
        final List<DiagnosticRow> rows = readDiagnosticRows(diagnosticFile);
        assertTrue("Expected dependent CTMC diagnostic rows", rows.size() > 0);

        boolean sawAvailableBeaglePreorder = false;
        for (DiagnosticRow row : rows) {
            final double exactDelta = parseDiagnosticDouble(row, "exactDelta");
            final double manualDelta = parseDiagnosticDouble(row, "manualDelta");
            assertFinite(exactDelta);
            assertFinite(manualDelta);
            assertEquals("manualDelta must match exactDelta on diagnostic row " + row.get("row"),
                    exactDelta, manualDelta, TOL);
            assertEquals("manualMinusExactDelta must be roundoff on diagnostic row " + row.get("row"),
                    0.0, parseDiagnosticDouble(row, "manualMinusExactDelta"), TOL);

            if ("available".equals(row.get("beaglePreorderStatus"))) {
                sawAvailableBeaglePreorder = true;
                final double beagleDelta = parseDiagnosticDouble(row, "beagleDelta");
                assertFinite(beagleDelta);
                assertEquals("beagleDelta must match exactDelta on diagnostic row " + row.get("row"),
                        exactDelta, beagleDelta, TOL);
                assertEquals("beagleMinusExactDelta must be roundoff on diagnostic row " + row.get("row"),
                        0.0, parseDiagnosticDouble(row, "beagleMinusExactDelta"), TOL);
            }
        }

        assertTrue("BEAGLE preorder diagnostics were not available; use the patched BEAGLE runtime for this test.",
                sawAvailableBeaglePreorder);
    }

    private static List<DiagnosticRow> readDiagnosticRows(final File diagnosticFile) throws Exception {
        final BufferedReader reader = new BufferedReader(new FileReader(diagnosticFile));
        try {
            final String header = reader.readLine();
            assertTrue("Diagnostic file must contain a header", header != null && header.length() > 0);

            final String[] columns = header.split("\t", -1);
            final List<DiagnosticRow> rows = new ArrayList<DiagnosticRow>();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.length() > 0) {
                    rows.add(new DiagnosticRow(columns, line.split("\t", -1)));
                }
            }
            return rows;
        } finally {
            reader.close();
        }
    }

    private static double parseDiagnosticDouble(final DiagnosticRow row, final String column) {
        return Double.parseDouble(row.get(column));
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

    private static void assertLogEquals(final double expected, final double actual) {
        assertLogEquals("Log values differ", expected, actual);
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static final class DiagnosticRow {
        final String[] header;
        final String[] values;

        private DiagnosticRow(final String[] header, final String[] values) {
            this.header = header;
            this.values = values;
        }

        private String get(final String column) {
            for (int i = 0; i < header.length; i++) {
                if (column.equals(header[i])) {
                    return values[i];
                }
            }
            throw new IllegalArgumentException("Missing diagnostic column " + column);
        }
    }

    private static final class Fixture {
        final TreeModel tree;
        final Parameter ctsRewards;
        final Parameter categories;
        final Parameter categoryCuts;
        final RewardsAwareCategoricalMixtureBranchRates branchRates;
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final TreeDataLikelihood independentLikelihood;
        final TreeDataLikelihood dependentLikelihood;

        private Fixture(final TreeModel tree,
                        final Parameter ctsRewards,
                        final Parameter categories,
                        final Parameter categoryCuts,
                        final RewardsAwareCategoricalMixtureBranchRates branchRates,
                        final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final TreeDataLikelihood independentLikelihood) {
            this(tree, ctsRewards, categories, categoryCuts, branchRates,
                    rewardsAwareBranchModel, independentLikelihood, null);
        }

        private Fixture(final TreeModel tree,
                        final Parameter ctsRewards,
                        final Parameter categories,
                        final Parameter categoryCuts,
                        final RewardsAwareCategoricalMixtureBranchRates branchRates,
                        final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final TreeDataLikelihood independentLikelihood,
                        final TreeDataLikelihood dependentLikelihood) {
            this.tree = tree;
            this.ctsRewards = ctsRewards;
            this.categories = categories;
            this.categoryCuts = categoryCuts;
            this.branchRates = branchRates;
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.independentLikelihood = independentLikelihood;
            this.dependentLikelihood = dependentLikelihood;
        }

        private Fixture withDependentLikelihood(final TreeDataLikelihood dependentLikelihood) {
            return new Fixture(tree, ctsRewards, categories, categoryCuts, branchRates,
                    rewardsAwareBranchModel, independentLikelihood, dependentLikelihood);
        }
    }
}
