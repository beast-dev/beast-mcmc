/*
 * BeagleRewardDependentCtmcEdgeEvidenceProviderTest.java
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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
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
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.BeagleRewardDependentCtmcEdgeEvidenceProvider;
import test.dr.math.MathTestCase;

/**
 * Regression tests for the BEAGLE preorder edge-evidence path used by the
 * reward-mixture indicator operator.
 *
 * @author Filippo Monti
 */
public class BeagleRewardDependentCtmcEdgeEvidenceProviderTest extends MathTestCase {

    private static final double TOL = 1.0e-9;

    public void testExternalEdgeEvidenceRatioMatchesFullLikelihoodRatio() {
        final Fixture fixture = createFixture();
        assertLocalEdgeEvidenceRatioMatchesFullLikelihoodRatio(fixture, fixture.tree.getExternalNode(0));
    }

    public void testInternalEdgeEvidenceRatioMatchesFullLikelihoodRatio() {
        final Fixture fixture = createFixture();
        assertLocalEdgeEvidenceRatioMatchesFullLikelihoodRatio(fixture, firstInternalNonRoot(fixture.tree));
    }

    private static void assertLocalEdgeEvidenceRatioMatchesFullLikelihoodRatio(final Fixture fixture,
                                                                               final NodeRef branch) {
        final BeagleRewardDependentCtmcEdgeEvidenceProvider provider =
                new BeagleRewardDependentCtmcEdgeEvidenceProvider(fixture.treeDataLikelihood);

        final int branchNodeNumber = branch.getNumber();
        final int branchParameterIndex = fixture.branchRateModel.getParameterIndexFromNode(branch);

        final double rawReward1 = 0.35;
        final double rawReward2 = 1.10;

        fixture.ctsRewards.setParameterValue(branchParameterIndex, rawReward1);
        fixture.treeDataLikelihood.makeDirty();
        final double logLikelihood1 = fixture.treeDataLikelihood.getLogLikelihood();

        provider.prepare();
        final double logEvidence1 = provider.logEvidence(branchNodeNumber, rawReward1);
        final double logEvidence2 = provider.logEvidence(branchNodeNumber, rawReward2);

        fixture.ctsRewards.setParameterValue(branchParameterIndex, rawReward2);
        fixture.treeDataLikelihood.makeDirty();
        final double logLikelihood2 = fixture.treeDataLikelihood.getLogLikelihood();

        assertFinite(logLikelihood1);
        assertFinite(logLikelihood2);
        assertFinite(logEvidence1);
        assertFinite(logEvidence2);
        assertEquals(logLikelihood1 - logLikelihood2, logEvidence1 - logEvidence2, TOL);
    }

    private static Fixture createFixture() {
        final TreeModel tree = createThreeTipTree();
        final SubstitutionModel substitutionModel = createNucleotideSubstitutionModel();
        final SitePatterns patterns = createSitePatterns();
        final GammaSiteRateModel siteRateModel = new GammaSiteRateModel("dependentRewardEvidenceRateModel");

        final Parameter ctsRewards = new Parameter.Default("dependentCtsRewards", new double[]{0.7, 0.9, 1.1, 1.3});
        final Parameter indicator = new Parameter.Default("dependentIndicator", new double[]{0.0, 0.0, 0.0, 0.0});
        final Parameter atomIndices = new Parameter.Default("dependentAtomIndices", new double[]{0.0, 1.0, 2.0, 3.0});
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("dependentRewardRates", new double[]{0.25, 0.75, 1.25, 1.75}),
                null,
                new Parameter.Default("dependentRewardRatesInternal", new double[0]),
                new Parameter.Default("dependentRewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
        final RewardsAwareMixtureBranchRates branchRateModel = new RewardsAwareMixtureBranchRates(
                tree,
                ctsRewards,
                indicator,
                atomIndices,
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );

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
        final TreeDataLikelihood treeDataLikelihood =
                new TreeDataLikelihood(delegate, tree, branchRateModel);

        return new Fixture(tree, treeDataLikelihood, branchRateModel, ctsRewards);
    }

    private static SitePatterns createSitePatterns() {
        final DataType dataType = Nucleotides.INSTANCE;
        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        alignment.addSequence(sequence("a", "ACGT", dataType));
        alignment.addSequence(sequence("b", "TGCA", dataType));
        alignment.addSequence(sequence("c", "CAGT", dataType));
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
        return new DefaultTreeModel("dependentRewardEvidenceTree", tree);
    }

    private static SubstitutionModel createNucleotideSubstitutionModel() {
        final DataType dataType = Nucleotides.INSTANCE;
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.25, 0.25, 0.25, 0.25});
        return new HKY(new Parameter.Default("dependentKappa", 2.0), frequencyModel);
    }

    private static NodeRef firstInternalNonRoot(final Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isExternal(node) && !tree.isRoot(node)) {
                return node;
            }
        }
        throw new IllegalArgumentException("Tree has no internal non-root branch");
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static final class Fixture {
        final TreeModel tree;
        final TreeDataLikelihood treeDataLikelihood;
        final RewardsAwareMixtureBranchRates branchRateModel;
        final Parameter ctsRewards;

        Fixture(final TreeModel tree,
                final TreeDataLikelihood treeDataLikelihood,
                final RewardsAwareMixtureBranchRates branchRateModel,
                final Parameter ctsRewards) {
            this.tree = tree;
            this.treeDataLikelihood = treeDataLikelihood;
            this.branchRateModel = branchRateModel;
            this.ctsRewards = ctsRewards;
        }
    }
}
