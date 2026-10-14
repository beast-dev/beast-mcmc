/*
 * CompleteHistoryLoggingTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.basta.BastaLikelihood;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.GenericBastaLikelihoodDelegate;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.markovjumps.MarkovJumpsLogger;
import dr.evomodel.treelikelihood.utilities.CompleteHistoryLogger;
import dr.evomodelxml.treedatalikelihood.markovjumps.MarkovJumpsParserUtils;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Stochastic-mapping traits and the loggers that read them: markovJumpsLogger and completeHistoryLogger for MASCOT,
 * and completeHistoryLogger for BASTA.
 */
public class CompleteHistoryLoggingTest extends TestCase {

    private static final String TAG = "location";

    public void testMascotTraitsDescribeOneDraw() throws Exception {
        MathUtils.setSeed(31);
        MascotMappingFixtures.Built built = MascotMappingFixtures.c3().build(0.02);
        MascotLikelihood likelihood = built.likelihood;
        Parameter fromZero = new Parameter.Default("fromZero", new double[]{0, 1, 1, 0, 0, 0, 0, 0, 0});
        Parameter timeInTwo = new Parameter.Default("timeInTwo", new double[]{0, 0, 1});
        likelihood.enableStochasticMapping(TAG, Collections.singletonList(fromZero),
                Collections.singletonList(timeInTwo));

        TreeModel tree = built.tree;
        TreeTrait states = likelihood.getTreeTrait(TAG);
        TreeTrait history = likelihood.getTreeTrait(TAG + MascotLikelihood.HISTORY_SUFFIX);
        TreeTrait count = likelihood.getTreeTrait(TAG + MascotLikelihood.COUNT_SUFFIX);
        TreeTrait raw = likelihood.getTreeTrait(TAG + MascotLikelihood.RAW_SUFFIX);
        TreeTrait jumpsFromZero = likelihood.getTreeTrait(TAG + ".fromZero");
        TreeTrait rewardInTwo = likelihood.getTreeTrait(TAG + ".timeInTwo");
        TreeTrait total = likelihood.getTreeTrait("allTransitions");
        assertNotNull(states);
        assertNotNull(history);
        assertNotNull(count);
        assertNotNull(jumpsFromZero);
        assertNotNull(rewardInTwo);
        assertNotNull(likelihood.getTreeTrait("history"));
        assertNotNull(total);
        assertFalse(raw.getLoggable());

        for (int step = 0; step < 20; step++) {
            nextLogStep(likelihood);
            double sum = 0.0;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                int state = ((int[]) states.getTrait(tree, node))[0];
                assertEquals(likelihood.getStochasticMap().getNodeState(node.getNumber()), state);
                assertEquals("\"" + state + "\"", states.getTraitString(tree, node));
                if (tree.isRoot(node)) {
                    assertNull(history.getTrait(tree, node));
                    continue;
                }
                StateHistory h = ((StateHistory[]) raw.getTrait(tree, node))[0];
                assertEquals(state, h.getStartingState());
                assertEquals(h.toStringChanges(1, built.dataType), ((String[]) history.getTrait(tree, node))[0]);
                assertEquals(h.getNumberOfJumps(), ((double[]) count.getTrait(tree, node))[0], 0.0);
                int[] jumps = h.getJumpCounts();
                assertEquals(jumps[1] + jumps[2], ((double[]) jumpsFromZero.getTrait(tree, node))[0], 0.0);
                assertEquals(h.getWaitingTimes()[2], ((double[]) rewardInTwo.getTrait(tree, node))[0], 0.0);
                sum += h.getNumberOfJumps();
            }
            assertEquals(sum, (Double) total.getTrait(tree, null), 0.0);
        }
    }

    /** Every logger in one log step reads the same draw; the next log step draws afresh. */
    public void testMascotLoggersAgreeWithinAStepAndRedrawBetweenSteps() throws Exception {
        MathUtils.setSeed(32);
        MascotMappingFixtures.Built built = MascotMappingFixtures.c3().build(0.02);
        MascotLikelihood likelihood = built.likelihood;
        likelihood.enableStochasticMapping(TAG, new ArrayList<>(), new ArrayList<>());

        CompleteHistoryLogger completeHistory = new CompleteHistoryLogger(likelihood, null, true, true);
        LogColumn[] historyColumns = completeHistory.getColumns();
        MarkovJumpsLogger jumps = new MarkovJumpsLogger(TAG + MascotLikelihood.COUNT_SUFFIX, likelihood,
                likelihood.getTreeTrait(TAG + MascotLikelihood.COUNT_SUFFIX), null, null, null, false);
        NumberColumn jumpColumn = (NumberColumn) jumps.getColumns()[0];

        Set<String> distinct = new HashSet<>();
        for (int step = 0; step < 30; step++) {
            nextLogStep(likelihood);
            String totalCount = historyColumns[0].getFormatted();
            String events = historyColumns[1].getFormatted();
            assertEquals(events, historyColumns[1].getFormatted());
            int logged = Integer.parseInt(events.substring(events.lastIndexOf(' ') + 1));
            assertEquals((double) logged, Double.parseDouble(totalCount), 0.0);
            assertEquals((double) logged, jumpColumn.getDoubleValue(), 0.0);
            assertEquals(likelihood.getStochasticMap().getTotalJumpCount(), logged);
            distinct.add(events);
        }
        assertTrue("draws did not change between log steps", distinct.size() > 5);
    }

    public void testMascotRejectsBadMappingRequests() throws Exception {
        MascotMappingFixtures.Built built = MascotMappingFixtures.c1().build(0.02);
        MascotLikelihood likelihood = built.likelihood;
        try {
            likelihood.enableStochasticMapping(TAG, Collections.singletonList(
                    new Parameter.Default("wrong", new double[]{1, 0, 0})), new ArrayList<>());
            fail("count register with dimension K accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            likelihood.enableStochasticMapping(TAG, Collections.singletonList(
                    new Parameter.Default("count", new double[]{0, 1, 1, 0})), new ArrayList<>());
            fail("register colliding with " + TAG + ".count accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            likelihood.enableStochasticMapping(MascotLikelihood.DEFAULT_MODE_STATE_TAG_NAME, new ArrayList<>(),
                    new ArrayList<>());
            fail("tag colliding with an existing trait accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        likelihood.enableStochasticMapping(TAG, new ArrayList<>(), new ArrayList<>());
        try {
            likelihood.enableStochasticMapping("other", new ArrayList<>(), new ArrayList<>());
            fail("second mapping accepted");
        } catch (IllegalStateException expected) {
            // expected
        }
    }

    /** BASTA's complete-history add-on now also feeds completeHistoryLogger. */
    public void testBastaCompleteHistoryLogger() throws Exception {
        MathUtils.setSeed(33);
        GeneralDataType dataType = new GeneralDataType(new String[]{"0", "1", "2"});
        TreeModel tree = new DefaultTreeModel("treeModel",
                new NewickImporter("(((A:1.0,B:1.0):1.2,(C:1.5,D:1.5):0.7):0.8,(E:2.0,F:2.0):1.0);").importTree(null));
        Patterns patterns = new Patterns(dataType, tree);
        int[] tipStates = {0, 1, 2, 0, 1, 2};
        int[] pattern = new int[tree.getTaxonCount()];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = tipStates[tree.getTaxonId(i).charAt(0) - 'A'];
        }
        patterns.addPattern(pattern);
        ComplexSubstitutionModel migration = new ComplexSubstitutionModel("location.model", dataType,
                new FrequencyModel(dataType, new double[]{1.0 / 3, 1.0 / 3, 1.0 / 3}),
                new Parameter.Default("rates", new double[]{0.6, 0.4, 0.5, 0.3, 0.7, 0.2}));
        migration.setNormalization(false);
        migration.setScaleRatesByFrequencies(false);
        Parameter popSizes = new Parameter.Default("popSizes", new double[]{1.5, 2.0, 1.2});
        BastaLikelihood likelihood = new BastaLikelihood("basta", tree, patterns, migration,
                new StrictClockBranchRates(new Parameter.Default("clock", 1.0)),
                new ConstantPopulationSizeModel("constant", popSizes, 3, -1), 1,
                new GenericBastaLikelihoodDelegate("basta", tree, 3, true), true, dataType, TAG, false);
        assertTrue(Double.isFinite(likelihood.getLogLikelihood()));
        assertSame(tree, likelihood.getTreeModel());
        assertEquals("\"2\"", likelihood.formattedState(new int[]{2}));

        MarkovJumpsParserUtils.addCompleteHistory(tree, likelihood, TAG);
        CompleteHistoryLogger completeHistory = new CompleteHistoryLogger(likelihood, null, true, true);
        LogColumn[] columns = completeHistory.getColumns();
        TreeTrait count = likelihood.getTreeTrait(TAG + ".count");
        List<String> seen = new ArrayList<>();
        for (int step = 0; step < 10; step++) {
            // BASTA redraws after a model change (or restore), not after makeDirty().
            popSizes.setParameterValue(0, 1.5 + 0.01 * step);
            String events = columns[1].getFormatted();
            int logged = Integer.parseInt(events.substring(events.lastIndexOf(' ') + 1));
            double sum = 0.0;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    sum += ((double[]) count.getTrait(tree, node))[0];
                }
            }
            assertEquals(sum, Double.parseDouble(columns[0].getFormatted()), 0.0);
            assertEquals(sum, logged, 0.0);
            seen.add(events);
        }
        assertTrue(new HashSet<>(seen).size() > 1);
    }

    /** What the MCMC does between log steps: store, propose nothing that touches MASCOT, accept. */
    private static void nextLogStep(MascotLikelihood likelihood) {
        likelihood.storeModelState();
        likelihood.acceptModelState();
    }
}
