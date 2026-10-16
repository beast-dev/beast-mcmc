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
