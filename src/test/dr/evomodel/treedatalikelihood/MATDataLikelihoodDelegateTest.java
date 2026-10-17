/*
 * MATDataLikelihoodDelegateTest.java
 *
 * Copyright © 2026, the BEAST Development Team.
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

package test.dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.parsimony.FitchParsimony;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.operators.FixedHeightSubtreePruneRegraftOperator;
import dr.evomodel.operators.SubtreeJumpOperator;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.MAT.DefaultMutationStore;
import dr.evomodel.treedatalikelihood.MAT.MATDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.Mutation;
import dr.evomodel.treedatalikelihood.MAT.MutationHistoryResampler;
import dr.evomodel.treedatalikelihood.MAT.MutationStore;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * Tests the MAT likelihood, its incremental parsimony and the mutation history resampler.
 */
public class MATDataLikelihoodDelegateTest extends TraceCorrelationAssert {

    private static final double CLOCK_RATE = 2.5;

    private SitePatterns patterns;
    private TreeModel treeModel;
    private HKY hky;
    private Parameter siteRates;
    private Parameter clockRate;
    private StrictClockBranchRates branchRates;
    private MATDataLikelihoodDelegate delegate;
    private TreeDataLikelihood likelihood;

    public MATDataLikelihoodDelegateTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(666);

        // primate sequences with some gaps, Ns and an ambiguity code to exercise ambiguities
        String[][] taxaSequence = new String[2][];
        taxaSequence[0] = PRIMATES_TAXON_SEQUENCE[0].clone();
        taxaSequence[1] = PRIMATES_TAXON_SEQUENCE[1].clone();
        taxaSequence[1][0] = mask(taxaSequence[1][0], 0, 40, '-');
        taxaSequence[1][1] = mask(taxaSequence[1][1], 20, 60, 'N');
        taxaSequence[1][2] = mask(taxaSequence[1][2], 30, 50, '?');
        taxaSequence[1][3] = mask(taxaSequence[1][3], 100, 101, 'R');
        createAlignment(taxaSequence, Nucleotides.INSTANCE);

        treeModel = createPrimateTreeModel();
        patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.CONSTANT_ONLY);

        FrequencyModel frequencyModel = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(new double[]{0.3, 0.2, 0.15, 0.35}));
        hky = new HKY(new Parameter.Default(4.0), frequencyModel);

        double[] rates = new double[patterns.getPatternCount()];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = 0.5 + (i % 7) * 0.25;
        }
        siteRates = new Parameter.Default(rates);
        clockRate = new Parameter.Default(CLOCK_RATE);
        branchRates = new StrictClockBranchRates(clockRate);

        delegate = new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky),
                null, branchRates, siteRates);
        likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
    }

    private static String mask(String sequence, int from, int to, char c) {
        StringBuilder sb = new StringBuilder(sequence);
        for (int i = from; i < to; i++) {
            sb.setCharAt(i, c);
        }
        return sb.toString();
    }

    public void testMutationCountMatchesFitch() {
        // the MAT treats partial ambiguities (R) as complete ambiguities so compare with Fitch on data with R -> N
        String[][] taxaSequence = new String[2][];
        taxaSequence[0] = PRIMATES_TAXON_SEQUENCE[0].clone();
        taxaSequence[1] = new String[taxaSequence[0].length];
        for (int i = 0; i < taxaSequence[0].length; i++) {
            taxaSequence[1][i] = alignment.getAlignedSequenceString(i).replace('R', 'N');
        }
        createAlignment(taxaSequence, Nucleotides.INSTANCE);
        SitePatterns unambiguousPatterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.CONSTANT_ONLY);

        assertEquals(new FitchParsimony(unambiguousPatterns, false).getScore(treeModel), weightedMutationCount(), 1E-10);
        assertTrue(delegate.getMutationStore().getAmbiguityCount() > 0);
    }

    /**
     * Each variable site needs its own pattern, so compressed variable patterns must be rejected, and the
     * constant-only compression must give exactly the same likelihood as uncompressed patterns.
     */
    public void testPatternCompression() {
        try {
            new MATDataLikelihoodDelegate(treeModel, new SitePatterns(alignment, null, 0, -1, 1, true),
                    new HomogeneousBranchModel(hky), null, branchRates, null);
            fail("compressed variable patterns should be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }

        double[] logL = new double[2];
        int i = 0;
        for (SitePatterns.CompressionType compression : new SitePatterns.CompressionType[]{
                SitePatterns.CompressionType.UNCOMPRESSED, SitePatterns.CompressionType.CONSTANT_ONLY}) {
            SitePatterns sitePatterns = new SitePatterns(alignment, null, 0, -1, 1, true, compression);
            MATDataLikelihoodDelegate d = new MATDataLikelihoodDelegate(treeModel, sitePatterns,
                    new HomogeneousBranchModel(hky), null, branchRates, null);
            logL[i++] = new TreeDataLikelihood(d, treeModel, branchRates).getLogLikelihood();
        }
        assertEquals(logL[0], logL[1], 1E-8);
    }

    public void testLikelihoodMatchesDirectPathDensity() {
        assertEquals(directLogLikelihood(), likelihood.getLogLikelihood(), 1E-6);

        // site rates and clock rate changes must propagate
        siteRates.setParameterValue(3, 4.0);
        clockRate.setParameterValue(0, 0.7);
        assertEquals(directLogLikelihood(), likelihood.getLogLikelihood(), 1E-6);
    }

    /**
     * After any sequence of topology moves (some rejected), incremental parsimony must give the
     * same MAT as building it from scratch on the current tree.
     */
    public void testIncrementalParsimonyMatchesRebuild() {
        List<AbstractTreeOperator> operators = new ArrayList<AbstractTreeOperator>();
        operators.add(new SubtreeLeapOperator(treeModel, 1.0, 0.05, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                AdaptationMode.ADAPTATION_OFF, 0.234));
        operators.add(new SubtreeJumpOperator(treeModel, 1.0, 1.0, 0.234, true, AdaptationMode.ADAPTATION_OFF));

        likelihood.getLogLikelihood();
        int changes = 0;
        for (int iteration = 0; iteration < 300; iteration++) {
            likelihood.storeModelState();
            String before = treeModel.getNewick();
            operators.get(iteration % operators.size()).doOperation();
            if (!treeModel.getNewick().equals(before)) {
                changes++;
            }
            double logL = likelihood.getLogLikelihood();

            MATDataLikelihoodDelegate fresh = new MATDataLikelihoodDelegate(treeModel, patterns,
                    new HomogeneousBranchModel(hky), null, branchRates, siteRates);
            assertSameStates(fresh.getMutationStore(), delegate.getMutationStore());
            assertEquals(new TreeDataLikelihood(fresh, treeModel, branchRates).getLogLikelihood(), logL, 1E-8);

            if (MathUtils.nextBoolean()) {
                likelihood.restoreModelState();
            } else {
                likelihood.acceptModelState();
            }
        }
        assertTrue(changes > 50);
        assertTrue(delegate.getRebuildCount() > 0);
    }

    /**
     * The resampler's proposal density for a block must be the full conditional of the block's
     * states, which is checked by enumerating every assignment for each pattern and scoring it with
     * the delegate's per-site likelihoods.
     */
    public void testResamplerMatchesEnumeratedConditional() {
        // a large block: parent of human plus the path from human's sibling (chimp+bonobo) to siamang
        checkResamplerAgainstEnumeration("human", "siamang", false);
        // a local block: parent of chimp plus the path from bonobo to human, which leaves variable sites
        // with no events around it
        checkResamplerAgainstEnumeration("chimp", "human", true);
    }

    /**
     * Without per-site rates the resampler filters in linear space, so check that path too.
     */
    public void testResamplerMatchesEnumeratedConditionalWithoutSiteRates() {
        delegate = new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky),
                null, branchRates, null);
        likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
        // a large block: parent of human plus the path from human's sibling (chimp+bonobo) to siamang
        checkResamplerAgainstEnumeration("human", "siamang", false);
        // a local block: parent of chimp plus the path from bonobo to human, which leaves variable sites
        // with no events around it
        checkResamplerAgainstEnumeration("chimp", "human", true);
    }

    /**
     * Checks the resampler's full conditional for the block of a move of the subtree at the given tip to
     * the branch above the given destination.
     */
    private void checkResamplerAgainstEnumeration(String moved, String destination, boolean expectNoEventSites) {
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate);
        DefaultMutationStore store = (DefaultMutationStore) delegate.getMutationStore();

        NodeRef node = findTip(moved);
        NodeRef parent = treeModel.getParent(node);
        int[] block = resampler.getBlock(parent, otherChild(parent, node), findTip(destination));
        assertTrue(block.length >= 2);

        // draw a fresh set of states so that we are not just testing the parsimony reconstruction
        likelihood.storeModelState();
        resampler.sampleBlock(block, true);
        updateBranches(store, block);
        likelihood.makeDirty();
        likelihood.getLogLikelihood();
        likelihood.acceptModelState();

        double resamplerLogProbability = resampler.sampleBlock(block, false);

        double enumeratedLogProbability = 0.0;
        int patternCount = patterns.getPatternCount();
        double[] original = likelihood.getDataLikelihoodDelegate().getSiteLogLikelihoods();
        byte[][] originalStates = new byte[block.length][];
        for (int i = 0; i < block.length; i++) {
            originalStates[i] = store.getStates(block[i]).clone();
        }

        likelihood.storeModelState();
        int[] assignment = new int[block.length];
        for (int p = 0; p < patternCount; p++) {
            if (delegate.isConstantPattern(p)) {
                // constant sites are assumed to have no mutations and are not part of the sampled state
                continue;
            }
            // nodes ambiguous at this pattern keep state -1; the others range over all states
            List<Integer> free = new ArrayList<Integer>();
            for (int i = 0; i < block.length; i++) {
                if (originalStates[i][p] >= 0) {
                    free.add(i);
                }
            }
            if (free.isEmpty()) {
                continue;
            }
            double logTotal = Double.NEGATIVE_INFINITY;
            int combinations = (int) Math.pow(4, free.size());
            for (int k = 0; k < combinations; k++) {
                int code = k;
                for (int i : free) {
                    assignment[i] = code % 4;
                    code /= 4;
                }
                for (int i : free) {
                    store.getWritableStates(block[i])[p] = (byte) assignment[i];
                }
                updateBranches(store, block);
                likelihood.makeDirty();
                likelihood.getLogLikelihood();
                double siteLogL = likelihood.getDataLikelihoodDelegate().getSiteLogLikelihoods()[p];
                logTotal = logAdd(logTotal, siteLogL);
            }
            for (int i : free) {
                store.getWritableStates(block[i])[p] = originalStates[i][p];
            }
            enumeratedLogProbability += original[p] - logTotal;
        }
        likelihood.restoreModelState();

        assertEquals(enumeratedLogProbability, resamplerLogProbability, 1E-8);

        // the fast path for patterns that are uniform around the block must score states identically
        byte[] classes = resampler.getPatternClasses(block);
        int[] classCounts = new int[4];
        for (byte c : classes) {
            classCounts[c]++;
        }
        String counts = java.util.Arrays.toString(classCounts);
        assertTrue(counts, classCounts[MutationHistoryResampler.HAS_SUBSTITUTIONS] > 0);
        if (expectNoEventSites) {
            assertTrue(counts, classCounts[MutationHistoryResampler.NO_EVENTS] > 0);
        }
        assertTrue(counts, classCounts[MutationHistoryResampler.CONSTANT] > 0);
        assertEquals(resamplerLogProbability, resampler.sampleBlock(block, false, classes), 1E-8);
    }

    /**
     * SPR moves with the resampler attached must leave a consistent MAT whose likelihood matches the
     * direct path density, and rejected moves must restore it exactly.
     */
    public void testResamplingMovesKeepMATConsistent() {
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate);
        List<AbstractTreeOperator> operators = new ArrayList<AbstractTreeOperator>();
        operators.add(new SubtreeLeapOperator(treeModel, 1.0, 0.05, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                AdaptationMode.ADAPTATION_OFF, 0.234));
        operators.add(new SubtreeJumpOperator(treeModel, 1.0, 1.0, 0.234, true, AdaptationMode.ADAPTATION_OFF));
        operators.add(new FixedHeightSubtreePruneRegraftOperator(treeModel, 1.0));
        for (AbstractTreeOperator operator : operators) {
            operator.setSubtreeMoveHook(resampler);
        }

        likelihood.getLogLikelihood();
        for (int iteration = 0; iteration < 300; iteration++) {
            double logLBefore = likelihood.getLogLikelihood();
            String statesBefore = describeStates();
            likelihood.storeModelState();

            double logHR = operators.get(iteration % operators.size()).doOperation();
            assertFalse(Double.isNaN(logHR));
            double logL = likelihood.getLogLikelihood();
            assertEquals(directLogLikelihood(), logL, 1E-6);

            if (MathUtils.nextBoolean()) {
                likelihood.restoreModelState();
                assertEquals(statesBefore, describeStates());
                assertEquals(logLBefore, likelihood.getLogLikelihood(), 1E-10);
            } else {
                likelihood.acceptModelState();
            }
        }
        assertEquals("the resampler should keep the MAT in step with the topology", 0, delegate.getRebuildCount());
    }

    /**
     * After resampling moves the Fitch sets are stale; a topology change by an operator without the
     * resampler must then rebuild the MAT exactly as parsimony from scratch would.
     */
    public void testParsimonyAfterResamplingMatchesRebuild() {
        SubtreeLeapOperator resampling = new SubtreeLeapOperator(treeModel, 1.0, 0.05,
                SubtreeLeapOperator.DistanceKernelType.NORMAL, AdaptationMode.ADAPTATION_OFF, 0.234);
        resampling.setSubtreeMoveHook(new MutationHistoryResampler(delegate));
        SubtreeJumpOperator plain = new SubtreeJumpOperator(treeModel, 1.0, 1.0, 0.234, true, AdaptationMode.ADAPTATION_OFF);

        likelihood.getLogLikelihood();
        for (int iteration = 0; iteration < 50; iteration++) {
            likelihood.storeModelState();
            resampling.doOperation();
            likelihood.getLogLikelihood();
            likelihood.acceptModelState();
        }
        assertEquals(0, delegate.getRebuildCount());

        likelihood.storeModelState();
        plain.doOperation();
        double logL = likelihood.getLogLikelihood();
        assertEquals(1, delegate.getRebuildCount());

        MATDataLikelihoodDelegate fresh = new MATDataLikelihoodDelegate(treeModel, patterns,
                new HomogeneousBranchModel(hky), null, branchRates, siteRates);
        assertSameStates(fresh.getMutationStore(), delegate.getMutationStore());
        assertEquals(new TreeDataLikelihood(fresh, treeModel, branchRates).getLogLikelihood(), logL, 1E-8);

        // and rejecting the move restores the sampled MAT
        likelihood.restoreModelState();
        assertEquals(directLogLikelihood(), likelihood.getLogLikelihood(), 1E-6);
    }

    /**
     * The hybrid resampler only resamples sites that vary around the block; for the move to be valid
     * the reverse move must resample exactly the same sites.
     */
    public void testResampledSitesAreSymmetric() {
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate);
        String[][] moves = {{"human", "gorilla"}, {"gorilla", "human"}, {"chimp", "orangutan"}, {"orangutan", "bonobo"}, {"siamang", "chimp"}};
        int checked = 0;
        for (String[] move : moves) {
            NodeRef x = findTip(move[0]);
            NodeRef destination = findTip(move[1]);
            NodeRef parent = treeModel.getParent(x);
            NodeRef sibling = otherChild(parent, x);
            if (treeModel.isRoot(parent) || destination == sibling || treeModel.getParent(destination) == parent) {
                continue;
            }
            double oldHeight = treeModel.getNodeHeight(parent);

            likelihood.storeModelState();
            resampler.beforeSubtreeMove(x, sibling, destination);
            regraft(x, destination, regraftHeight(destination, treeModel.getNodeHeight(x)));
            resampler.afterSubtreeMove();
            int[] forward = resampler.getLastResampledSites();
            likelihood.getLogLikelihood();

            // and back again
            resampler.beforeSubtreeMove(x, destination, sibling);
            regraft(x, sibling, oldHeight);
            resampler.afterSubtreeMove();
            int[] reverse = resampler.getLastResampledSites();
            likelihood.getLogLikelihood();
            likelihood.restoreModelState();

            java.util.Arrays.sort(forward);
            java.util.Arrays.sort(reverse);
            assertTrue(move[0] + " -> " + move[1] + ": " + java.util.Arrays.toString(forward) + " vs " + java.util.Arrays.toString(reverse),
                    java.util.Arrays.equals(forward, reverse));
            assertTrue(forward.length > 0);
            checked++;
        }
        assertTrue(checked >= 3);
    }

    /**
     * With site selection ALL and no per-site rates, the resampler handles the sites with no events
     * around the block together (drawing the few that change by thinning). Its Hastings ratio must still
     * be that of the exact Gibbs block update: the full conditional of the old states under the old
     * topology over that of the new states under the new one, computed here over every pattern.
     */
    public void testSparseResamplingMatchesDenseConditional() {
        delegate = new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky),
                null, branchRates, null);
        likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
        // the direct path density (directLogLikelihood) then uses a relative rate of 1 for every site
        siteRates = new Parameter.Default(patterns.getPatternCount(), 1.0);
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate, MutationHistoryResampler.SiteSelection.ALL);

        int moves = 0;
        int movesChangingUniformSites = 0;
        for (int iteration = 0; iteration < 2000 && moves < 300; iteration++) {
            NodeRef x = treeModel.getNode(MathUtils.nextInt(treeModel.getNodeCount()));
            NodeRef destination = treeModel.getNode(MathUtils.nextInt(treeModel.getNodeCount()));
            if (treeModel.isRoot(x) || treeModel.isRoot(treeModel.getParent(x)) || treeModel.isRoot(destination)) {
                continue;
            }
            NodeRef parent = treeModel.getParent(x);
            NodeRef sibling = otherChild(parent, x);
            if (destination == parent || destination == sibling || isAncestorOrSelf(x, destination) ||
                    treeModel.getNodeHeight(treeModel.getParent(destination)) <= treeModel.getNodeHeight(x)) {
                // (the regrafted parent must fit between x and the destination's parent)
                continue;
            }

            likelihood.getLogLikelihood();
            likelihood.storeModelState();
            int[] block = resampler.getBlock(parent, sibling, destination);
            double denseOld = resampler.sampleBlock(block, false);

            resampler.beforeSubtreeMove(x, sibling, destination);
            regraft(x, destination, regraftHeight(destination, treeModel.getNodeHeight(x)));
            double logRatio = resampler.afterSubtreeMove();

            double denseNew = resampler.sampleBlock(block, false);
            assertEquals("move " + moves, denseOld - denseNew, logRatio, 1E-8);
            assertEquals(directLogLikelihood(), likelihood.getLogLikelihood(), 1E-6);
            if (resampler.getLastChangedUniformSiteCount() > 0) {
                movesChangingUniformSites++;
            }

            if (MathUtils.nextBoolean()) {
                likelihood.acceptModelState();
            } else {
                likelihood.restoreModelState();
            }
            moves++;
        }
        assertEquals(300, moves);
        // the sites with no events around the block must sometimes change, or this tests nothing new
        assertTrue("moves changing uniform sites: " + movesChangingUniformSites, movesChangingUniformSites > 0);
    }

    /**
     * The sites with no events around the block, and those with only ambiguities there, are drawn by
     * thinning; check that each changes with its exact probability (from the full conditional of its block
     * states, enumerated with the dense scorer) by repeating the same move many times.
     */
    public void testSparseResamplingDrawsUniformSitesExactly() {
        // a high rate, so that sites with no events around the block change often enough to count
        clockRate.setParameterValue(0, 40.0);
        delegate = new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky),
                null, branchRates, null);
        likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate, MutationHistoryResampler.SiteSelection.ALL);
        DefaultMutationStore store = (DefaultMutationStore) delegate.getMutationStore();

        NodeRef x = findTip("chimp");
        NodeRef destination = findTip("human");
        NodeRef parent = treeModel.getParent(x);
        NodeRef sibling = otherChild(parent, x);
        double height = regraftHeight(destination, treeModel.getNodeHeight(x));
        // (restoring the tree replaces its node objects, so the nodes are looked up by number each time)
        final int xNumber = x.getNumber();
        final int siblingNumber = sibling.getNumber();
        final int destinationNumber = destination.getNumber();
        likelihood.getLogLikelihood();
        int[] block = resampler.getBlock(parent, sibling, destination);
        byte[] classes = resampler.getPatternClasses(block);
        List<Integer> uniformSites = new ArrayList<Integer>();
        int[] uniformState = new int[patterns.getPatternCount()];
        int[] neighbourhood = new int[]{x.getNumber(), sibling.getNumber(), destination.getNumber(),
                treeModel.getParent(destination).getNumber(), parent.getNumber()};
        int noEventSites = 0;
        int ambiguitySites = 0;
        for (int p = 0; p < patterns.getPatternCount(); p++) {
            if (classes[p] != MutationHistoryResampler.NO_EVENTS && classes[p] != MutationHistoryResampler.AMBIGUITIES_ONLY) {
                continue;
            }
            // every unambiguous node around the block has the same state
            int a = -1;
            for (int n : neighbourhood) {
                a = Math.max(a, store.getState(n, p));
            }
            for (int n : block) {
                a = Math.max(a, store.getState(n, p));
            }
            if (a >= 0) {
                uniformSites.add(p);
                uniformState[p] = a;
                if (classes[p] == MutationHistoryResampler.NO_EVENTS) {
                    noEventSites++;
                } else {
                    ambiguitySites++;
                }
            }
        }
        assertTrue("sites with no events: " + noEventSites, noEventSites > 10);
        assertTrue("sites with only ambiguities: " + ambiguitySites, ambiguitySites > 10);

        // the exact probability that each such site's block states change, under the new topology
        likelihood.storeModelState();
        resampler.beforeSubtreeMove(x, sibling, destination);
        regraft(x, destination, height);
        resampler.afterSubtreeMove();
        double[] probabilityChanged = new double[patterns.getPatternCount()];
        boolean[][] free = new boolean[patterns.getPatternCount()][];
        double expectedChanges = 0.0;
        for (int p : uniformSites) {
            // the block nodes with data below them at this site (in the new topology); the others stay ambiguous
            free[p] = new boolean[block.length];
            int freeCount = 0;
            for (int i = 0; i < block.length; i++) {
                free[p][i] = store.getState(block[i], p) >= 0;
                if (free[p][i]) {
                    freeCount++;
                }
            }
            double logUnchanged = Double.NaN;
            double logTotal = Double.NEGATIVE_INFINITY;
            int combinations = (int) Math.pow(4, freeCount);
            for (int k = 0; k < combinations; k++) {
                int code = k;
                boolean unchanged = true;
                for (int i = 0; i < block.length; i++) {
                    if (!free[p][i]) {
                        continue;
                    }
                    int state = code % 4;
                    code /= 4;
                    store.getWritableStates(block[i])[p] = (byte) state;
                    unchanged &= state == uniformState[p];
                }
                double logQ = resampler.sampleBlock(block, false);
                logTotal = logAdd(logTotal, logQ);
                if (unchanged) {
                    logUnchanged = logQ;
                }
            }
            for (int i = 0; i < block.length; i++) {
                if (free[p][i]) {
                    store.getWritableStates(block[i])[p] = (byte) uniformState[p];
                }
            }
            probabilityChanged[p] = -Math.expm1(logUnchanged - logTotal);
            expectedChanges += probabilityChanged[p];
        }
        likelihood.restoreModelState();

        // repeat the move and count the changes
        final int repeats = 2000;
        int[] changes = new int[patterns.getPatternCount()];
        int totalChanges = 0;
        for (int r = 0; r < repeats; r++) {
            x = treeModel.getNode(xNumber);
            sibling = treeModel.getNode(siblingNumber);
            destination = treeModel.getNode(destinationNumber);
            likelihood.storeModelState();
            resampler.beforeSubtreeMove(x, sibling, destination);
            regraft(x, destination, height);
            resampler.afterSubtreeMove();
            for (int p : uniformSites) {
                for (int i = 0; i < block.length; i++) {
                    int state = store.getState(block[i], p);
                    assertEquals("ambiguity of block node " + i + " at site " + p, free[p][i], state >= 0);
                    if (free[p][i] && state != uniformState[p]) {
                        changes[p]++;
                        totalChanges++;
                        break;
                    }
                }
            }
            likelihood.restoreModelState();
        }

        // the total number of changes is a sum of independent Bernoullis
        double mean = repeats * expectedChanges;
        double variance = 0.0;
        for (int p : uniformSites) {
            variance += repeats * probabilityChanged[p] * (1.0 - probabilityChanged[p]);
        }
        assertTrue("expected " + mean + " changes, observed " + totalChanges, mean > 50);
        assertEquals("changes of the sites with no events around the block", mean, totalChanges, 4 * Math.sqrt(variance));
        // and site by site
        for (int p : uniformSites) {
            double siteMean = repeats * probabilityChanged[p];
            double siteSd = Math.sqrt(siteMean * (1.0 - probabilityChanged[p]));
            assertEquals("site " + p, siteMean, changes[p], 5 * siteSd + 1);
        }
    }

    /**
     * Data with regions in which only one or two taxa have data, so that moves have sites at which the
     * moved subtree is the only data along part of its path (hot sites), with chains of different lengths.
     */
    private void setUpSparseData() {
        String[][] taxaSequence = new String[2][];
        taxaSequence[0] = PRIMATES_TAXON_SEQUENCE[0].clone();
        taxaSequence[1] = PRIMATES_TAXON_SEQUENCE[1].clone();
        for (int t = 0; t < taxaSequence[1].length; t++) {
            String sequence = taxaSequence[1][t];
            if (t != 4) {
                sequence = mask(sequence, 200, 260, 'N'); // orangutan only
            }
            if (t != 4 && t != 5) {
                sequence = mask(sequence, 300, 340, 'N'); // orangutan and siamang
            }
            if (t != 0 && t != 1) {
                sequence = mask(sequence, 400, 430, 'N'); // human and chimp
            }
            if (t != 3) {
                sequence = mask(sequence, 500, 520, '-'); // gorilla only
            }
            if (t == 2 || t == 3) {
                sequence = mask(sequence, 600, 640, 'N'); // all but bonobo and gorilla
            }
            taxaSequence[1][t] = sequence;
        }
        createAlignment(taxaSequence, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();
        patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.CONSTANT_ONLY);
        branchRates = new StrictClockBranchRates(clockRate);
        delegate = new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky),
                null, branchRates, null);
        likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
        // the direct path density (directLogLikelihood) then uses a relative rate of 1 for every site
        siteRates = new Parameter.Default(patterns.getPatternCount(), 1.0);
    }

    /**
     * The local block: the log Hastings ratio of random moves equals log q(old) - log q(new), where q is the
     * full conditional of the block's states (the regrafted parent's, and those of every node whose
     * presence the move changes) given all other states, found by enumerating the block at each site.
     * States outside the block must not change, and the MAT must stay consistent.
     */
    public void testLocalBlockMatchesEnumeratedConditional() {
        setUpSparseData();
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate,
                MutationHistoryResampler.SiteSelection.ALL, true);

        int moves = 0;
        int chainSites = 0;
        int longChainSites = 0;
        for (int iteration = 0; iteration < 5000 && moves < 200; iteration++) {
            NodeRef x = treeModel.getNode(MathUtils.nextInt(treeModel.getNodeCount()));
            NodeRef destination = treeModel.getNode(MathUtils.nextInt(treeModel.getNodeCount()));
            if (treeModel.isRoot(x) || treeModel.isRoot(treeModel.getParent(x)) || treeModel.isRoot(destination)) {
                continue;
            }
            NodeRef parent = treeModel.getParent(x);
            NodeRef sibling = otherChild(parent, x);
            if (destination == parent || destination == sibling || isAncestorOrSelf(x, destination) ||
                    treeModel.getNodeHeight(treeModel.getParent(destination)) <= treeModel.getNodeHeight(x)) {
                continue;
            }

            likelihood.getLogLikelihood();
            likelihood.storeModelState();
            String statesBefore = describeStates();
            double logLBefore = likelihood.getLogLikelihood();

            int p = parent.getNumber();
            int[] oldParents = currentParents();
            int[] newParents = oldParents.clone();
            newParents[sibling.getNumber()] = oldParents[p];
            newParents[p] = oldParents[destination.getNumber()];
            newParents[destination.getNumber()] = p;
            boolean[][] oldPresence = presence(oldParents);
            boolean[][] newPresence = presence(newParents);
            int[][] oldStates = snapshotStates();
            double logQOld = 0.0;
            for (int site : variableSites()) {
                int[] block = localBlock(p, site, oldPresence, newPresence);
                logQOld += logConditional(site, block, oldPresence, oldStates, oldParents);
            }

            resampler.beforeSubtreeMove(x, sibling, destination);
            regraft(x, destination, regraftHeight(destination, treeModel.getNodeHeight(x)));
            double logRatio = resampler.afterSubtreeMove();

            assertTrue(java.util.Arrays.equals(newParents, currentParents()));
            int[][] newStates = snapshotStates();
            double logQNew = 0.0;
            for (int site : variableSites()) {
                int[] block = localBlock(p, site, oldPresence, newPresence);
                boolean[] inBlock = new boolean[treeModel.getNodeCount()];
                for (int n : block) {
                    inBlock[n] = true;
                }
                for (int n = 0; n < treeModel.getNodeCount(); n++) {
                    if (!inBlock[n]) {
                        assertEquals("node " + n + " outside the block at site " + site, oldStates[n][site], newStates[n][site]);
                    }
                }
                if (block.length > 1) {
                    chainSites++;
                    if (block.length > 2) {
                        longChainSites++;
                    }
                }
                logQNew += logConditional(site, block, newPresence, newStates, newParents);
            }
            assertEquals("move " + moves, logQOld - logQNew, logRatio, 1E-6);
            assertEquals(directLogLikelihood(), likelihood.getLogLikelihood(), 1E-6);

            if (MathUtils.nextBoolean()) {
                likelihood.acceptModelState();
            } else {
                likelihood.restoreModelState();
                assertEquals(statesBefore, describeStates());
                assertEquals(logLBefore, likelihood.getLogLikelihood(), 1E-10);
            }
            moves++;
        }
        assertEquals(200, moves);
        assertTrue("sites with chains: " + chainSites, chainSites > 0);
        assertTrue("sites with chains longer than one node: " + longChainSites, longChainSites > 0);
    }

    /**
     * The local block's states are drawn from their full conditional: repeat one move (with long chains)
     * many times and compare the frequency of each site's drawn block states with the enumerated conditional.
     */
    public void testLocalBlockDrawsExactly() {
        clockRate.setParameterValue(0, 40.0);
        setUpSparseData();
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate,
                MutationHistoryResampler.SiteSelection.ALL, true);

        final int xNumber = findTip("orangutan").getNumber();
        final int destinationNumber = findTip("human").getNumber();
        likelihood.getLogLikelihood();
        NodeRef x = treeModel.getNode(xNumber);
        NodeRef parent = treeModel.getParent(x);
        final int p = parent.getNumber();
        final int siblingNumber = otherChild(parent, x).getNumber();
        int[] oldParents = currentParents();
        int[] newParents = oldParents.clone();
        newParents[siblingNumber] = oldParents[p];
        newParents[p] = oldParents[destinationNumber];
        newParents[destinationNumber] = p;
        boolean[][] oldPresence = presence(oldParents);
        boolean[][] newPresence = presence(newParents);
        final int[] sites = variableSites();
        final double height = regraftHeight(treeModel.getNode(destinationNumber), treeModel.getNodeHeight(x));

        final int repeats = 3000;
        java.util.Map<String, Integer> observed = new java.util.HashMap<String, Integer>();
        java.util.Map<String, Double> expected = new java.util.HashMap<String, Double>();
        int longChainSites = 0;
        for (int r = 0; r < repeats; r++) {
            likelihood.getLogLikelihood();
            likelihood.storeModelState();
            x = treeModel.getNode(xNumber);
            resampler.beforeSubtreeMove(x, treeModel.getNode(siblingNumber), treeModel.getNode(destinationNumber));
            regraft(x, treeModel.getNode(destinationNumber), height);
            resampler.afterSubtreeMove();
            int[][] states = snapshotStates();
            for (int site : sites) {
                int[] block = present(localBlock(p, site, oldPresence, newPresence), newPresence, site);
                if (r == 0) {
                    if (block.length > 2) {
                        longChainSites++;
                    }
                    // the exact distribution of this site's block states (the other states are the same each time)
                    double logZ = logConditionalNormaliser(site, block, states, newParents);
                    int[] assignment = states(site, block, states).clone();
                    for (int code = 0; code < (1 << (2 * block.length)); code++) {
                        for (int i = 0; i < block.length; i++) {
                            states[block[i]][site] = (code >> (2 * i)) & 3;
                        }
                        double probability = Math.exp(siteLogLikelihood(site, states, newParents) - logZ);
                        expected.put(site + ":" + java.util.Arrays.toString(states(site, block, states)), probability);
                    }
                    for (int i = 0; i < block.length; i++) {
                        states[block[i]][site] = assignment[i];
                    }
                }
                String key = site + ":" + java.util.Arrays.toString(states(site, block, states));
                Integer count = observed.get(key);
                observed.put(key, count == null ? 1 : count + 1);
            }
            likelihood.restoreModelState();
        }
        assertTrue("sites with chains longer than one node: " + longChainSites, longChainSites > 0);

        int checked = 0;
        for (java.util.Map.Entry<String, Double> e : expected.entrySet()) {
            double mean = repeats * e.getValue();
            Integer count = observed.get(e.getKey());
            int n = count == null ? 0 : count;
            assertEquals(e.getKey(), mean, n, 5 * Math.sqrt(mean * (1 - e.getValue())) + 1);
            if (mean > 5 && mean < repeats - 5) {
                checked++;
            }
        }
        for (String key : observed.keySet()) {
            assertTrue("drawn states with no probability: " + key, expected.containsKey(key));
        }
        assertTrue("block states with intermediate probabilities: " + checked, checked > 20);
    }

    /**
     * The local block with the subtree moves themselves (including moves at the root): the MAT stays
     * consistent and restores exactly.
     */
    public void testLocalBlockMovesKeepMATConsistent() {
        setUpSparseData();
        MutationHistoryResampler resampler = new MutationHistoryResampler(delegate,
                MutationHistoryResampler.SiteSelection.ALL, true);
        List<AbstractTreeOperator> operators = new ArrayList<AbstractTreeOperator>();
        operators.add(new SubtreeLeapOperator(treeModel, 1.0, 0.05, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                AdaptationMode.ADAPTATION_OFF, 0.234));
        operators.add(new SubtreeJumpOperator(treeModel, 1.0, 1.0, 0.234, true, AdaptationMode.ADAPTATION_OFF));
        operators.add(new FixedHeightSubtreePruneRegraftOperator(treeModel, 1.0));
        for (AbstractTreeOperator operator : operators) {
            operator.setSubtreeMoveHook(resampler);
        }

        likelihood.getLogLikelihood();
        for (int iteration = 0; iteration < 600; iteration++) {
            double logLBefore = likelihood.getLogLikelihood();
            String statesBefore = describeStates();
            likelihood.storeModelState();

            double logHR = operators.get(iteration % operators.size()).doOperation();
            assertFalse(Double.isNaN(logHR));
            double logL = likelihood.getLogLikelihood();
            assertEquals(directLogLikelihood(), logL, 1E-6);

            if (MathUtils.nextBoolean()) {
                likelihood.restoreModelState();
                assertEquals(statesBefore, describeStates());
                assertEquals(logLBefore, likelihood.getLogLikelihood(), 1E-10);
            } else {
                likelihood.acceptModelState();
            }
        }
        assertEquals("the resampler should keep the MAT in step with the topology", 0, delegate.getRebuildCount());
    }

    private int[] currentParents() {
        int[] parents = new int[treeModel.getNodeCount()];
        for (int n = 0; n < parents.length; n++) {
            NodeRef parent = treeModel.getParent(treeModel.getNode(n));
            parents[n] = parent == null ? -1 : parent.getNumber();
        }
        return parents;
    }

    private int[] variableSites() {
        List<Integer> sites = new ArrayList<Integer>();
        for (int site = 0; site < patterns.getPatternCount(); site++) {
            if (!delegate.isConstantPattern(site)) {
                sites.add(site);
            }
        }
        int[] result = new int[sites.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = sites.get(i);
        }
        return result;
    }

    /**
     * @return [node][site]: whether a tip below the node (in the topology given by parents) has data
     */
    private boolean[][] presence(int[] parents) {
        MutationStore store = delegate.getMutationStore();
        boolean[][] present = new boolean[parents.length][patterns.getPatternCount()];
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            int tip = treeModel.getExternalNode(i).getNumber();
            for (int site = 0; site < patterns.getPatternCount(); site++) {
                if (store.getState(tip, site) >= 0) {
                    for (int n = tip; n >= 0 && !present[n][site]; n = parents[n]) {
                        present[n][site] = true;
                    }
                }
            }
        }
        return present;
    }

    private int[][] snapshotStates() {
        MutationStore store = delegate.getMutationStore();
        int[][] states = new int[treeModel.getNodeCount()][patterns.getPatternCount()];
        for (int n = 0; n < states.length; n++) {
            for (int site = 0; site < states[n].length; site++) {
                states[n][site] = store.getState(n, site);
            }
        }
        return states;
    }

    /**
     * @return the local block at a site: the regrafted parent and every node whose presence the move changes
     */
    private static int[] localBlock(int p, int site, boolean[][] oldPresence, boolean[][] newPresence) {
        List<Integer> block = new ArrayList<Integer>();
        block.add(p);
        for (int n = 0; n < oldPresence.length; n++) {
            if (n != p && oldPresence[n][site] != newPresence[n][site]) {
                block.add(n);
            }
        }
        int[] result = new int[block.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = block.get(i);
        }
        return result;
    }

    private static int[] present(int[] nodes, boolean[][] presence, int site) {
        List<Integer> result = new ArrayList<Integer>();
        for (int n : nodes) {
            if (presence[n][site]) {
                result.add(n);
            }
        }
        int[] array = new int[result.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = result.get(i);
        }
        return array;
    }

    private static int[] states(int site, int[] nodes, int[][] states) {
        int[] result = new int[nodes.length];
        for (int i = 0; i < nodes.length; i++) {
            result[i] = states[nodes[i]][site];
        }
        return result;
    }

    /**
     * @return the log full conditional of the block's states at a site (given all other states), in the
     * topology given by parents: the site's likelihood less its log sum over the states of the block's
     * nodes that have data
     */
    private double logConditional(int site, int[] block, boolean[][] presence, int[][] states, int[] parents) {
        for (int n : block) {
            if (!presence[n][site]) {
                assertEquals("node " + n + " without data at site " + site, -1, states[n][site]);
            }
        }
        int[] free = present(block, presence, site);
        return siteLogLikelihood(site, states, parents) - logConditionalNormaliser(site, free, states, parents);
    }

    private double logConditionalNormaliser(int site, int[] free, int[][] states, int[] parents) {
        int[] saved = states(site, free, states);
        double logZ = Double.NEGATIVE_INFINITY;
        for (int code = 0; code < (1 << (2 * free.length)); code++) {
            for (int i = 0; i < free.length; i++) {
                states[free[i]][site] = (code >> (2 * i)) & 3;
            }
            logZ = logAdd(logZ, siteLogLikelihood(site, states, parents));
        }
        for (int i = 0; i < free.length; i++) {
            states[free[i]][site] = saved[i];
        }
        return logZ;
    }

    /**
     * @return the log likelihood of one site's states (single hits, a relative rate of 1) in the topology
     * given by parents, with the current node heights
     */
    private double siteLogLikelihood(int site, int[][] states, int[] parents) {
        double[] q = new double[16];
        hky.getInfinitesimalMatrix(q);
        double[] pi = hky.getFrequencyModel().getFrequencies();
        double rate = clockRate.getParameterValue(0);
        double logL = 0.0;
        for (int n = 0; n < parents.length; n++) {
            int b = states[n][site];
            if (parents[n] < 0) {
                if (b >= 0) {
                    logL += Math.log(pi[b]);
                }
                continue;
            }
            if (b < 0) {
                continue;
            }
            int a = states[parents[n]][site];
            assertTrue("unambiguous child below ambiguous parent", a >= 0);
            double length = (treeModel.getNodeHeight(treeModel.getNode(parents[n])) -
                    treeModel.getNodeHeight(treeModel.getNode(n))) * rate;
            logL += singleHitLogDensity(q, a, b, length, 1.0);
        }
        return logL;
    }

    private static double singleHitLogDensity(double[] q, int a, int b, double length, double nu) {
        if (a == b) {
            return -nu * -q[a * 4 + a] * length;
        }
        double qa = -q[a * 4 + a], qb = -q[b * 4 + b];
        int steps = 2000;
        double h = length / steps, integral = 0;
        for (int k = 0; k <= steps; k++) {
            double s = k * h;
            double f = nu * q[a * 4 + b] * Math.exp(-nu * qa * s - nu * qb * (length - s));
            integral += f * (k == 0 || k == steps ? 1 : (k % 2 == 1 ? 4 : 2));
        }
        return Math.log(integral * h / 3);
    }

    private boolean isAncestorOrSelf(NodeRef ancestor, NodeRef node) {
        for (NodeRef n = node; n != null; n = treeModel.getParent(n)) {
            if (n == ancestor) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return a height for the regrafted parent on the branch above the destination (and above x)
     */
    private double regraftHeight(NodeRef destination, double minimum) {
        double lower = Math.max(treeModel.getNodeHeight(destination), minimum);
        NodeRef parent = treeModel.getParent(destination);
        double upper = parent == null ? lower + 0.01 : treeModel.getNodeHeight(parent);
        return 0.5 * (lower + upper);
    }

    /**
     * Prunes x and its parent and regrafts them on the branch above the destination at the given height.
     */
    private void regraft(NodeRef x, NodeRef destination, double height) {
        NodeRef parent = treeModel.getParent(x);
        NodeRef sibling = otherChild(parent, x);
        NodeRef grandParent = treeModel.getParent(parent);
        NodeRef destinationParent = treeModel.getParent(destination);
        treeModel.beginTreeEdit();
        treeModel.removeChild(parent, sibling);
        treeModel.removeChild(grandParent, parent);
        treeModel.addChild(grandParent, sibling);
        treeModel.removeChild(destinationParent, destination);
        treeModel.addChild(parent, destination);
        treeModel.addChild(destinationParent, parent);
        treeModel.setNodeHeight(parent, height);
        treeModel.endTreeEdit();
    }

    private double weightedMutationCount() {
        MutationStore store = delegate.getMutationStore();
        double count = 0;
        for (int n = 0; n < store.getNodeCount(); n++) {
            for (Mutation m : store.getMutationsForNode(n)) {
                if (!m.isAmbiguity()) {
                    count += patterns.getPatternWeight(m.getPattern());
                }
            }
        }
        return count;
    }

    private void updateBranches(DefaultMutationStore store, int[] block) {
        for (int n : block) {
            store.updateBranch(n);
            NodeRef node = treeModel.getNode(n);
            for (int i = 0; i < treeModel.getChildCount(node); i++) {
                store.updateBranch(treeModel.getChild(node, i).getNumber());
            }
        }
    }

    private static double logAdd(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) {
            return b;
        }
        double max = Math.max(a, b);
        return max + Math.log(Math.exp(a - max) + Math.exp(b - max));
    }

    private String describeStates() {
        MutationStore store = delegate.getMutationStore();
        StringBuilder sb = new StringBuilder(treeModel.getNewick());
        for (int n = 0; n < store.getNodeCount(); n++) {
            sb.append('|').append(store.getParent(n)).append(':');
            for (int p = 0; p < store.getPatternCount(); p++) {
                sb.append(store.getState(n, p));
            }
            sb.append(':').append(store.getMutationsForNode(n).length);
        }
        return sb.toString();
    }

    /**
     * The store's node summaries must equal those computed from the states.
     */
    private void assertCompositionsMatchStates() {
        DefaultMutationStore store = (DefaultMutationStore) delegate.getMutationStore();
        for (int n = 0; n < treeModel.getNodeCount(); n++) {
            double[] expected = new double[4];
            for (int p = 0; p < patterns.getPatternCount(); p++) {
                int a = store.getState(n, p);
                if (a >= 0) {
                    expected[a] += patterns.getPatternWeight(p);
                }
            }
            double[] actual = store.getComposition(n);
            for (int a = 0; a < 4; a++) {
                assertEquals("composition of node " + n + " state " + a, expected[a], actual[a], 1E-9);
            }
        }
    }

    private void assertSameStates(MutationStore expected, MutationStore actual) {
        for (int n = 0; n < expected.getNodeCount(); n++) {
            assertEquals(expected.getParent(n), actual.getParent(n));
            for (int p = 0; p < expected.getPatternCount(); p++) {
                assertEquals("node " + n + " pattern " + p, expected.getState(n, p), actual.getState(n, p));
            }
            assertEquals(expected.getMutationsForNode(n).length, actual.getMutationsForNode(n).length);
        }
    }

    private NodeRef findTip(String id) {
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef tip = treeModel.getExternalNode(i);
            if (treeModel.getNodeTaxon(tip).getId().equals(id)) {
                return tip;
            }
        }
        throw new IllegalArgumentException(id);
    }

    private NodeRef otherChild(NodeRef parent, NodeRef child) {
        NodeRef c = treeModel.getChild(parent, 0);
        return c == child ? treeModel.getChild(parent, 1) : c;
    }

    /**
     * Sums, per site and branch, the log of the CTMC path density between the stored node states
     * (with any mutation time integrated by Simpson's rule), also checking that the store's states
     * are consistent with its topology, its branch lists and the data.
     */
    private double directLogLikelihood() {
        MutationStore store = delegate.getMutationStore();
        int nodeCount = treeModel.getNodeCount();
        int patternCount = patterns.getPatternCount();

        double[] q = new double[16];
        hky.getInfinitesimalMatrix(q);
        double[] pi = hky.getFrequencyModel().getFrequencies();

        NodeRef root = treeModel.getRoot();
        assertEquals(root.getNumber(), store.getRoot());
        assertCompositionsMatchStates();

        double logL = 0;
        for (int p = 0; p < patternCount; p++) {
            int a = store.getState(root.getNumber(), p);
            if (a >= 0) {
                logL += patterns.getPatternWeight(p) * Math.log(pi[a]);
            }
        }

        for (int n = 0; n < nodeCount; n++) {
            NodeRef node = treeModel.getNode(n);
            NodeRef parent = treeModel.getParent(node);
            if (parent == null) {
                continue;
            }
            assertEquals(parent.getNumber(), store.getParent(n));
            double length = (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node))
                    * branchRates.getBranchRate(treeModel, node);
            int mutationCount = 0;
            boolean childPresent = treeModel.isExternal(node);
            for (int p = 0; p < patternCount; p++) {
                int a = store.getState(parent.getNumber(), p);
                int b = store.getState(n, p);
                double nu = siteRates.getParameterValue(p);
                double siteLogL;
                if (b < 0) {
                    siteLogL = 0;
                } else {
                    assertTrue("unambiguous child below ambiguous parent", a >= 0);
                    if (a == b) {
                        siteLogL = -nu * -q[a * 4 + a] * length;
                    } else {
                        mutationCount++;
                        double qa = -q[a * 4 + a], qb = -q[b * 4 + b];
                        int steps = 2000;
                        double h = length / steps, integral = 0;
                        for (int k = 0; k <= steps; k++) {
                            double s = k * h;
                            double f = nu * q[a * 4 + b] * Math.exp(-nu * qa * s - nu * qb * (length - s));
                            integral += f * (k == 0 || k == steps ? 1 : (k % 2 == 1 ? 4 : 2));
                        }
                        siteLogL = Math.log(integral * h / 3);
                    }
                }
                logL += patterns.getPatternWeight(p) * siteLogL;
            }
            int storedMutationCount = 0;
            for (Mutation m : store.getMutationsForNode(n)) {
                if (!m.isAmbiguity()) {
                    storedMutationCount++;
                }
            }
            assertEquals("mutation list of node " + n, mutationCount, storedMutationCount);
            // the events of a branch are in pattern order (substitutions, then ambiguities)
            dr.evomodel.treedatalikelihood.MAT.BranchMutations events = store.getBranchMutations(n);
            for (int k = 1; k < events.size(); k++) {
                if (k != events.getSubstitutionCount()) {
                    assertTrue("events of node " + n + " out of pattern order", events.getPattern(k - 1) <= events.getPattern(k));
                }
            }

            // an internal node is ambiguous exactly where all its children are
            if (!childPresent) {
                for (int p = 0; p < patternCount; p++) {
                    boolean anyChild = false;
                    for (int i = 0; i < treeModel.getChildCount(node); i++) {
                        anyChild |= store.getState(treeModel.getChild(node, i).getNumber(), p) >= 0;
                    }
                    assertEquals("node " + n + " pattern " + p + " state " + store.getState(n, p) + " children " + treeModel.getChildCount(node), anyChild, store.getState(n, p) >= 0);
                }
            }
        }

        // tips must agree with the data where unambiguous and be ambiguous (-1) otherwise, including R
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef tip = treeModel.getExternalNode(i);
            int taxonIndex = patterns.getTaxonIndex(treeModel.getNodeTaxon(tip).getId());
            for (int p = 0; p < patternCount; p++) {
                int observed = patterns.getPatternState(taxonIndex, p);
                assertEquals(observed < 4 ? observed : -1, store.getState(tip.getNumber(), p));
            }
        }
        return logL;
    }
}
