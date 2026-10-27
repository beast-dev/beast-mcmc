/*
 * TimeTreeSearchTest.java
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

package test.dr.app.treelime;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.PartitionSubstitutionModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.BranchRateWarmStart;
import dr.evomodel.tree.search.NniNeighbourhood;
import dr.evomodel.tree.search.RootedNniMove;
import dr.evomodel.tree.search.SearchLog;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.evomodel.tree.search.TopologyHash;
import dr.evomodel.tree.search.TrialEvaluator;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * The topology search's first version (plan 7 M1a) on XML TreeLime writes for dated Dengue sequences, with the
 * round maximizer and the search added by hand: V1 (the NNI neighbourhood), V2 (a trial leaves no trace), V3
 * (accepting by replay), V6 (the rate warm start), V7 (the objective never decreases), V8 (the constraint filter),
 * V9 (gradients after topology changes), V13 (determinism), plus outputs and parser refusals. Run with assertions on
 * (TreeDataLikelihood.storeState and TreeIntervals check invariants with assert).
 */
public class TimeTreeSearchTest extends TestCase {

    private static final String DENGUE = "examples/Data/Dengue4.env.nex";

    public TimeTreeSearchTest(String name) {
        super(name);
    }

    // V1: 2 (n - 2) rooted NNIs, each feasible after the move (repaired where needed), with and without tip dates
    public void testNeighbourhood() throws Exception {
        for (boolean dated : new boolean[]{true, false}) {
            Map<String, XMLObject> objects = parse(new Setup().dated(dated).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, null);
            assertEquals(2 * (tree.getExternalNodeCount() - 2), moves.size());
            int repaired = 0;
            Model model = joint(objects).getModel();
            for (RootedNniMove move : moves) {
                joint(objects).getLogLikelihood();
                repaired += move.needsRepair(tree) ? 1 : 0;
                model.storeModelState();
                move.apply(tree, null);
                assertValidHeights(tree);
                assertFalse(Double.isInfinite(joint(objects).getLogLikelihood()));
                model.restoreModelState();
            }
            assertTrue("some moves need a repair", repaired > 0);
        }
    }

    // V2: every trial (T0, T1, T2), for every move, leaves the tree, the parameters, the density and the gradients as
    // they were; a full recomputation agrees (verifyEvery)
    public void testTrialLeavesNoTrace() throws Exception {
        for (boolean relaxed : new boolean[]{false, true}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(relaxed).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            TimeTreeSearch search = search(objects);
            perturb(tree, 2);
            TrialEvaluator evaluator = new TrialEvaluator(tree, joint, roundMaximizer(objects),
                    new BranchRateWarmStart(tree, relaxed ? rates(objects) : null,
                            BranchRateWarmStart.Mode.PRESERVE_SUBSTITUTIONS),
                    new dr.evomodel.tree.search.LocalOptimizer(tree, relaxed ? rates(objects) : null, 2), 5, 1);
            GradientWrtParameterProvider gradient =
                    (GradientWrtParameterProvider) objects.get("gradient.nodeHeights").getNativeObject();

            final String newick = TreeUtils.newick(tree);
            final double f = joint.getLogLikelihood();
            final double[] parameters = roundParameters(evaluator);
            final double[] heightGradient = gradient.getGradientLogDensity().clone();
            final Likelihood coalescent = (Likelihood) objects.get("coalescent").getNativeObject();
            final double coalescentValue = coalescent.getLogLikelihood();

            int k = 0;
            for (RootedNniMove move : NniNeighbourhood.enumerate(tree, null)) {
                TrialEvaluator.Trial trial = evaluator.trial(move, Double.NEGATIVE_INFINITY, k++ % 5 == 0);
                assertFalse(Double.isNaN(trial.f1));
                assertEquals(newick, TreeUtils.newick(tree));
                assertEquals(f, joint.getLogLikelihood(), 1E-10 * Math.abs(f));
                assertArrayEquals(parameters, roundParameters(evaluator), 0.0);
                assertEquals(coalescentValue, coalescent.getLogLikelihood(), 1E-10 * Math.abs(coalescentValue));
                // a direct gradient call right after the restore, with no transform() first
                assertArrayEquals(heightGradient, gradient.getGradientLogDensity(), 1E-8);
            }
            assertTrue(evaluator.getVerifications() >= k);
            assertTrue(evaluator.getT2Count() > 0);
            assertNotNull(search);
        }
    }

    // V3: an accepted move, replayed with the recorded values, gives the trial's F; a full recomputation agrees
    public void testReplay() throws Exception {
        for (boolean relaxed : new boolean[]{false, true}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(relaxed).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            TrialEvaluator evaluator = new TrialEvaluator(tree, joint, roundMaximizer(objects),
                    new BranchRateWarmStart(tree, relaxed ? rates(objects) : null,
                            BranchRateWarmStart.Mode.PRESERVE_SUBSTITUTIONS),
                    new dr.evomodel.tree.search.LocalOptimizer(tree, relaxed ? rates(objects) : null, 2), 5, 0);
            for (int k = 0; k < 4; k++) {
                List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, null);
                TrialEvaluator.Trial trial = evaluator.trial(moves.get((5 * k + 1) % moves.size()),
                        Double.NEGATIVE_INFINITY, k % 2 == 1);
                final String before = TopologyHash.toHex(TopologyHash.of(tree));
                final double f = evaluator.accept(trial);
                assertEquals(trial.best, f, 1E-8 * Math.abs(f));
                assertFalse(before.equals(TopologyHash.toHex(TopologyHash.of(tree))));
                joint.makeDirty();
                assertEquals(f, joint.getLogLikelihood(), 1E-8 * Math.abs(f));
            }
        }
    }

    // V6: the warm start keeps rate times time on every rewritten branch, b and P included when P is repaired
    public void testWarmStartPreservesSubstitutions() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        ArbitraryBranchRates rates = rates(objects);
        Parameter raw = rates.getRateParameter();
        for (int i = 0; i < raw.getDimension(); i++) {
            raw.setParameterValue(i, 0.001 * (1.0 + 0.1 * (i % 7)));
        }
        BranchRateWarmStart warmStart = new BranchRateWarmStart(tree, rates,
                BranchRateWarmStart.Mode.PRESERVE_SUBSTITUTIONS);
        Model model = joint(objects).getModel();
        int repaired = 0;
        int plain = 0;
        for (RootedNniMove move : NniNeighbourhood.enumerate(tree, null)) {
            joint(objects).getLogLikelihood();
            final boolean repair = move.needsRepair(tree);
            final int b = otherChild(tree, move.getP(), move.getA());
            final int[] nodes = repair ? new int[]{move.getA(), move.getU(), b, move.getP()}
                    : new int[]{move.getA(), move.getU()};
            double[] before = new double[nodes.length];
            double[] oldRates = new double[nodes.length];
            for (int i = 0; i < nodes.length; i++) {
                NodeRef node = tree.getNode(nodes[i]);
                oldRates[i] = rates.getBranchRate(tree, node);
                before[i] = oldRates[i] * branchTime(tree, node);
            }
            model.storeModelState();
            move.apply(tree, warmStart);
            for (int i = 0; i < nodes.length; i++) {
                NodeRef node = tree.getNode(nodes[i]);
                final double rate = rates.getBranchRate(tree, node);
                final double factor = rate / oldRates[i];
                if (factor > 0.2 * (1 + 1E-9) && factor < 5.0 * (1 - 1E-9)) {
                    assertEquals(move + " node " + nodes[i], before[i], rate * branchTime(tree, node),
                            1E-9 * before[i]);
                } else {
                    // clamped
                    assertTrue(Math.abs(factor - 0.2) < 1E-9 || Math.abs(factor - 5.0) < 1E-9);
                }
            }
            model.restoreModelState();
            if (repair) {
                repaired++;
            } else {
                plain++;
            }
        }
        assertTrue(repaired > 0 && plain > 0);
    }

    // V7: from perturbed starts, under both clocks and with a gamma shape held at its bound, F never decreases
    // (every round's F after is at least its F before, the best only rises) and the search improves on the start
    public void testMonotone() throws Exception {
        for (Setup setup : new Setup[]{new Setup(), new Setup().relaxed(true), new Setup().gamma(true)}) {
            Map<String, XMLObject> objects = parse(setup.xml());
            if (objects == null) {
                return;
            }
            perturb(tree(objects), 3);
            TimeTreeSearch search = search(objects);
            search.run();
            double best = Double.NEGATIVE_INFINITY;
            double previous = Double.NEGATIVE_INFINITY;
            for (String row : search.getLog().getRows()) {
                String[] fields = row.split("\t", -1);
                double before = Double.parseDouble(fields[7]);
                double after = Double.parseDouble(fields[8]);
                double fBest = Double.parseDouble(fields[9]);
                assertTrue(row, after >= before - 1E-6 * Math.abs(before));
                assertTrue(row, fBest >= best - 1E-9 * Math.abs(fBest));
                assertTrue(row, after >= previous - 1E-6 * Math.abs(after));
                best = fBest;
                previous = after;
            }
            assertTrue(search.getAcceptedMoves() > 0);
            assertTrue(search.getFinalValue() >= search.getBestValue() - 1E-6 * Math.abs(search.getBestValue()));
            assertEquals(TimeTreeSearch.StopReason.CONVERGED, search.getStopReason());
            assertEquals(search.getFinalValue(), joint(objects).getLogLikelihood(), 1E-8 * Math.abs(search.getFinalValue()));
        }
    }

    // V8: the constraint filter refuses exactly the NNIs after which the monophyly statistic fails; the search keeps
    // the clade
    public void testConstraintFilter() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().clade(new int[]{1, 2, 3}).xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        Likelihood joint = joint(objects);
        TopologyConstraints constraints = new TopologyConstraints(tree, joint, null);
        assertEquals(1, constraints.size());
        assertTrue(constraints.areSatisfied(tree));
        Model model = joint.getModel();
        int refused = 0;
        for (RootedNniMove move : NniNeighbourhood.enumerate(tree, null)) {
            joint.getLogLikelihood();
            final boolean violates = constraints.violates(tree, move);
            model.storeModelState();
            move.apply(tree, null);
            assertEquals(move.describe(), violates, joint.getLogLikelihood() == Double.NEGATIVE_INFINITY);
            model.restoreModelState();
            refused += violates ? 1 : 0;
        }
        assertTrue(refused > 0);
        assertEquals(NniNeighbourhood.enumerate(tree, null).size() - refused,
                NniNeighbourhood.enumerate(tree, constraints).size());

        perturb(tree, 2, constraints);
        TimeTreeSearch search = search(objects);
        search.run();
        assertEquals(1, search.getConstraints().size());
        assertTrue(constraints.areSatisfied(tree));
        assertTrue(search.getFinalValue() > Double.NEGATIVE_INFINITY);
    }

    // V9: the gradients still match finite differences in the re-fits after topology changes (the round
    // maximizer's list is its latest call's). With the HMC relaxed clock finite differences are poor where rates
    // near 0 (as in TreeLimeMaximizationTest), so there only no gradient may be NaN
    public void testGradientsAfterMoves() throws Exception {
        for (boolean relaxed : new boolean[]{false, true}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(relaxed).checkGradients(true).xml());
            if (objects == null) {
                return;
            }
            perturb(tree(objects), 3);
            TimeTreeSearch search = search(objects);
            search.run();
            assertTrue(search.getAcceptedMoves() > 0);
            for (String id : new String[]{"roundMaximizer", "maximizer"}) {
                List<Double> errors = ((MaximizerWrtParameter) objects.get(id).getNativeObject()).getGradientCheckErrors();
                assertFalse(id, errors.isEmpty());
                int large = 0;
                for (double error : errors) {
                    assertFalse(Double.isNaN(error));
                    if (!(error < 1E-4)) {
                        large++;
                    }
                }
                if (!relaxed) {
                    // a missing or double-counted term fails nearly every check
                    List<Double> sorted = new ArrayList<Double>(errors);
                    Collections.sort(sorted);
                    assertTrue(id + ": median error " + sorted.get(sorted.size() / 2),
                            sorted.get(sorted.size() / 2) < 1E-5);
                    assertTrue(id + ": " + large + " of " + errors.size() + " checks above 1E-4",
                            large <= Math.max(1, errors.size() / 20));
                }
            }
        }
    }

    // V13: the same XML and start give the same log, row for row (but the seconds), and the same tree
    public void testDeterminism() throws Exception {
        String xml = new Setup().relaxed(true).search("seed=\"17\"").xml();
        if (xml == null) {
            return;
        }
        List<String> logs = new ArrayList<String>();
        List<String> trees = new ArrayList<String>();
        for (int run = 0; run < 2; run++) {
            Map<String, XMLObject> objects = parse(xml);
            perturb(tree(objects), 3);
            TimeTreeSearch search = search(objects);
            search.run();
            StringBuilder sb = new StringBuilder();
            for (String row : search.getLog().getRows()) {
                sb.append(row, 0, row.lastIndexOf('\t')).append('\n');
            }
            logs.add(sb.toString());
            trees.add(TreeUtils.newick(tree(objects)));
        }
        assertEquals(logs.get(0), logs.get(1));
        assertEquals(trees.get(0), trees.get(1));
    }

    // the log and the NEXUS tree files
    public void testOutputs() throws Exception {
        File directory = Files.createTempDirectory("timeTreeSearch").toFile();
        File treeFile = new File(directory, "run.search.trees");
        File logFile = new File(directory, "run.search.log");
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).search("treeFileName=\"" +
                treeFile.getAbsolutePath() + "\" logFileName=\"" + logFile.getAbsolutePath() + "\"").xml());
        if (objects == null) {
            return;
        }
        perturb(tree(objects), 2);
        TimeTreeSearch search = search(objects);
        search.run();

        List<String> lines = Files.readAllLines(logFile.toPath());
        assertEquals(String.join("\t", SearchLog.COLUMNS), lines.get(0));
        assertEquals(search.getLog().getRows().size() + 1, lines.size());

        Tree written = new NexusImporter(new FileReader(treeFile)).importNextTree();
        assertEquals(tree(objects).getExternalNodeCount(), written.getExternalNodeCount());
        assertEquals(TopologyHash.of(tree(objects)), TopologyHash.of(written));
        assertEquals(tree(objects).getNodeHeight(tree(objects).getRoot()),
                (Double) written.getNodeAttribute(written.getRoot(), "height"), 1E-9);
        String text = new String(Files.readAllBytes(treeFile.toPath()));
        assertTrue(text.contains("&F=" + search.getFinalValue()));
        assertTrue(text.contains("rate="));
        String report = search.getReport();
        assertTrue(report.contains("StopReason: converged"));
        assertTrue(report.contains("AcceptedMoves: " + search.getAcceptedMoves()));
    }

    // the parser's refusals
    public void testParserRefusals() throws Exception {
        String xml = new Setup().relaxed(true).xml();
        if (xml == null) {
            return;
        }
        assertRefused(xml.replace("<maximizeWrtParameter id=\"roundMaximizer\" maximizeOnParse=\"false\" " +
                        "nIterations=\"200\" printToScreen=\"false\" startAtCurrentState=\"true\" includeJacobian=\"false\"",
                "<maximizeWrtParameter id=\"roundMaximizer\" maximizeOnParse=\"false\" " +
                        "nIterations=\"200\" printToScreen=\"false\" startAtCurrentState=\"true\" includeJacobian=\"true\""),
                "Jacobian");
        assertRefused(xml.replace("<timeTreeSearch id=\"search\" searchOnParse=\"false\" printToScreen=\"false\" >\n" +
                        "  <treeModel idref=\"treeModel\"/>\n  <objective><joint idref=\"joint\"/></objective>",
                "<timeTreeSearch id=\"search\" searchOnParse=\"false\" printToScreen=\"false\" >\n" +
                        "  <treeModel idref=\"treeModel\"/>\n  <objective><likelihood idref=\"likelihood\"/></objective>"),
                "another density");
        assertRefused(xml.replace("<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"true\""), "already run");
        assertRefused(xml.replace("  <jointGradient idref=\"gradient.nodeHeights\"/>\n", ""), "node heights");
    }

    // -------------------------------------------------------------------------------------------------------------

    private static void assertRefused(String xml, String message) {
        try {
            parse(xml);
            fail("not refused: " + message);
        } catch (Exception e) {
            String text = e.toString();
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                text += " / " + cause;
            }
            assertTrue(text, text.contains(message));
        }
    }

    // the TreeLime XML plus the search blocks
    static final class Setup {
        boolean relaxed = false;
        boolean dated = true;
        boolean gamma = false;
        boolean checkGradients = false;
        int[] clade = null;
        String searchAttributes = "";

        Setup relaxed(boolean relaxed) {
            this.relaxed = relaxed;
            return this;
        }

        Setup dated(boolean dated) {
            this.dated = dated;
            return this;
        }

        Setup gamma(boolean gamma) {
            this.gamma = gamma;
            return this;
        }

        Setup checkGradients(boolean check) {
            this.checkGradients = check;
            return this;
        }

        Setup clade(int[] taxa) {
            this.clade = taxa;
            return this;
        }

        Setup search(String attributes) {
            this.searchAttributes = attributes;
            return this;
        }

        String xml() throws Exception {
            TreeLimeOptions options = TreeLimeHarness.importData(DENGUE);
            if (options == null) {
                return null;
            }
            if (dated) {
                TreeLimeHarness.setTipDatesFromNames(options);
            }
            if (relaxed) {
                for (PartitionClockModel clock : options.getPartitionClockModels()) {
                    clock.setClockType(ClockType.HMC_CLOCK);
                }
            }
            if (gamma) {
                for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                    model.setGammaHetero(true);
                }
            }
            if (clade != null) {
                dr.evolution.util.Taxa taxa = new dr.evolution.util.Taxa("clade");
                for (int i : clade) {
                    taxa.addTaxon(options.taxonList.getTaxon(i));
                }
                options.taxonSets.add(taxa);
                options.taxonSetsMono.put(taxa, true);
                options.taxonSetsIncludeStem.put(taxa, false);
                options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
            }
            options.maximizationOptions.printToScreen = false;
            String xml = TreeLimeHarness.generate(options);
            xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                    "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"" +
                            (checkGradients ? " checkGradientEvery=\"1\"" : ""));
            xml = xml.substring(0, xml.indexOf("<!-- The optimum"));
            final String rates = relaxed ? "branchRatesHMC.rates" : (dated ? "clock.rate" : null);
            StringBuilder sb = new StringBuilder(xml);
            if (rates != null) {
                assertTrue(xml.contains("<jointGradient id=\"gradient." + rates + "\">"));
                sb.append("<compoundParameter id=\"searchRates\"><parameter idref=\"").append(rates)
                        .append("\"/></compoundParameter>\n");
            }
            sb.append("<compoundGradient id=\"searchGradient\">\n");
            if (rates != null) {
                sb.append("  <jointGradient idref=\"gradient.").append(rates).append("\"/>\n");
            }
            sb.append("  <jointGradient idref=\"gradient.nodeHeights\"/>\n</compoundGradient>\n");
            sb.append("<maximizeWrtParameter id=\"roundMaximizer\" maximizeOnParse=\"false\" nIterations=\"200\" " +
                    "printToScreen=\"false\" startAtCurrentState=\"true\" includeJacobian=\"false\" lineSearchRestarts=\"5\"")
                    .append(checkGradients ? " checkGradientEvery=\"1\"" : "").append(">\n")
                    .append("  <objective><joint idref=\"joint\"/></objective>\n")
                    .append("  <compoundGradient idref=\"searchGradient\"/>\n")
                    .append("  <multivariateCompoundTransform>\n");
            if (rates != null) {
                sb.append("    <compoundTransform><transform type=\"log\"><parameter idref=\"").append(rates)
                        .append("\"/></transform><compoundParameter idref=\"searchRates\"/></compoundTransform>\n");
            }
            sb.append("    <nodeHeightTransform idref=\"nodeHeightTransform\"/>\n")
                    .append("  </multivariateCompoundTransform>\n</maximizeWrtParameter>\n");
            sb.append("<timeTreeSearch id=\"search\" searchOnParse=\"false\" printToScreen=\"false\" ")
                    .append(searchAttributes).append(">\n")
                    .append("  <treeModel idref=\"treeModel\"/>\n")
                    .append("  <objective><joint idref=\"joint\"/></objective>\n")
                    .append("  <roundMaximizer><maximizeWrtParameter idref=\"roundMaximizer\"/></roundMaximizer>\n")
                    .append("  <finalMaximizer><maximizeWrtParameter idref=\"maximizer\"/></finalMaximizer>\n")
                    .append("</timeTreeSearch>\n</beast>\n");
            return sb.toString();
        }
    }

    static Map<String, XMLObject> parse(String xml) throws Exception {
        if (xml == null) {
            return null;
        }
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    static TreeModel tree(Map<String, XMLObject> objects) {
        return (TreeModel) objects.get("treeModel").getNativeObject();
    }

    static Likelihood joint(Map<String, XMLObject> objects) {
        return (Likelihood) objects.get("joint").getNativeObject();
    }

    static TimeTreeSearch search(Map<String, XMLObject> objects) {
        return (TimeTreeSearch) objects.get("search").getNativeObject();
    }

    static MaximizerWrtParameter roundMaximizer(Map<String, XMLObject> objects) {
        return (MaximizerWrtParameter) objects.get("roundMaximizer").getNativeObject();
    }

    static ArbitraryBranchRates rates(Map<String, XMLObject> objects) {
        return (ArbitraryBranchRates) objects.get("branchRates").getNativeObject();
    }

    // k repaired rooted NNIs, at fixed positions in the move order
    static void perturb(TreeModel tree, int k) {
        perturb(tree, k, null);
    }

    static void perturb(TreeModel tree, int k, TopologyConstraints constraints) {
        for (int i = 0; i < k; i++) {
            List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, constraints);
            moves.get((7 * i + 3) % moves.size()).apply(tree, null);
        }
    }

    private static double[] roundParameters(TrialEvaluator evaluator) {
        List<Double> values = new ArrayList<Double>();
        for (Parameter parameter : evaluator.getRoundParameters()) {
            for (double value : parameter.getParameterValues()) {
                values.add(value);
            }
        }
        double[] array = new double[values.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = values.get(i);
        }
        return array;
    }

    private static int otherChild(TreeModel tree, int parent, int child) {
        NodeRef p = tree.getNode(parent);
        NodeRef first = tree.getChild(p, 0);
        return first.getNumber() == child ? tree.getChild(p, 1).getNumber() : first.getNumber();
    }

    private static double branchTime(TreeModel tree, NodeRef node) {
        return tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
    }

    private static void assertValidHeights(TreeModel tree) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            for (int c = 0; c < tree.getChildCount(node); c++) {
                assertTrue(tree.getNodeHeight(node) > tree.getNodeHeight(tree.getChild(node, c)));
            }
        }
    }

    private static void assertArrayEquals(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("element " + i, expected[i], actual[i], tolerance * Math.max(1.0, Math.abs(expected[i])));
        }
    }
}
