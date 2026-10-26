/*
 * TreeLimeMaximizationTest.java
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
import dr.app.treelime.options.MaximizationOptions;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.PartitionSubstitutionModel;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.app.treelime.types.StartingTreeType;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * The XML that TreeLime's generator writes for the Maximization tab: it parses, the maximizer optimizes the chosen
 * density over the chosen parameters, and it stops by converging. With analytic gradients every run also checks, at
 * every evaluation, that the gradient blocks add up to the derivative of the objective (a missing or double-counted
 * term shows up there), and the optimum must match the all-numeric run's.
 */
public class TreeLimeMaximizationTest extends TestCase {

    public TreeLimeMaximizationTest(String name) {
        super(name);
    }

    // the gradient checks go along random directions: fix them, so a run does not depend on what ran before
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        dr.math.MathUtils.setSeed(666);
    }

    // BEAUti's default model: HKY, strict clock, constant-size coalescent
    public void testMapOfTheDefaultModel() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        String xml = TreeLimeHarness.generate(quiet(options));

        assertTrue(xml.contains("<maximizeWrtParameter id=\"maximizer\" nIterations=\"200\""));
        assertTrue(xml.contains("includeJacobian=\"false\""));
        assertTrue(between(xml, "<objective>", "</objective>").contains("<joint idref=\"joint\"/>"));
        // kappa: the analytic substitution-model gradient plus the numeric prior term; popSize: the prior term only
        assertTrue(between(xml, "<jointGradient id=\"gradient.kappa\">", "</jointGradient>")
                .contains("<branchSubstitutionParameterGradient homogeneous=\"true\">"));
        assertTrue(xml.contains("<numericalGradient id=\"gradient.constant.popSize\">"));
        assertTrue(xml.contains("usePreOrder=\"true\" branchInfinitesimalDerivative=\"true\""));
        // kappa and the population size on the log scale; the frequencies (a simplex) are held fixed
        assertTrue(xml.contains("<transform type=\"log\">"));
        assertFalse(between(xml, "<compoundParameter id=\"optimizedParameters\">", "</compoundParameter>")
                .contains("frequencies"));
        assertFalse(xml.contains("<mcmc"));
        assertFalse(xml.contains("<operators"));

        Run run = run(xml);
        assertEquals("joint", run.maximizer.getLikelihood().getId());
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue("the joint density increased", run.end > run.start);
        assertTrue(run.parameter("kappa").getParameterValue(0) != run.startValue("kappa"));
        assertGradientChecked(run);

        // the all-numeric run reaches the same optimum; it holds the node heights and the frequencies fixed, so hold
        // them here too
        options.maximizationOptions.excluded.add("treeModel.rootHeight");
        options.maximizationOptions.excluded.add("frequencies");
        run = run(TreeLimeHarness.generate(options));
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        options.maximizationOptions.setGradientMode(MaximizationOptions.GradientMode.NUMERIC);
        String numericXml = TreeLimeHarness.generate(options);
        assertTrue(numericXml.contains("<densityWrtParameter>"));
        assertFalse(numericXml.contains("branchInfinitesimalDerivative"));
        Run numeric = run(numericXml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, numeric.maximizer.getReason());
        for (String id : new String[]{"kappa", "constant.popSize"}) {
            double analytic = run.parameter(id).getParameterValue(0);
            assertEquals(id, analytic, numeric.parameter(id).getParameterValue(0), 1E-4 * analytic);
        }
        assertEquals(run.end, numeric.end, 1E-6 * Math.abs(run.end));
    }

    public void testGammaShape() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setGammaHetero(true);
        }
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<jointGradient id=\"gradient.alpha\">", "</jointGradient>")
                .contains("<siteRateModelGradient>"));
        // the gamma shape's default lower bound
        assertTrue(xml.contains("<holdAtLowerBound value=\"0.05\">\n\t\t\t<parameter idref=\"alpha\"/>"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);

        // a blank bound writes none
        for (dr.app.treelime.options.Parameter parameter : options.selectParameters()) {
            if (parameter.getBaseName().equals("alpha")) {
                options.maximizationOptions.setLowerBound(parameter, null);
            }
        }
        assertFalse(TreeLimeHarness.generate(options).contains("<holdAtLowerBound"));
    }

    // a monophyletic taxon set that plain UPGMA splits: the UPGMA start is built with the set as a constraint
    public void testMonophyleticSetInTheStartingTree() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        dr.evolution.util.Taxa clade = new dr.evolution.util.Taxa("clade");
        clade.addTaxon(options.taxonList.getTaxon(0));
        clade.addTaxon(options.taxonList.getTaxon(7));
        options.taxonSets.add(clade);
        options.taxonSetsMono.put(clade, true);
        options.taxonSetsIncludeStem.put(clade, false);
        options.taxonSetsTreeModel.put(clade, options.getPartitionTreeModels().get(0));
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<upgmaTree>", "</upgmaTree>").contains("<constraint>"));

        java.util.Set<String> pair = new java.util.HashSet<String>(java.util.Arrays.asList(
                clade.getTaxon(0).getId(), clade.getTaxon(1).getId()));
        // without the constraint UPGMA splits the pair
        String unconstrained = xml.replaceAll("(?s)<constraint>.*?</constraint>", "");
        Run plain = run(unconstrained);
        assertFalse(dr.evolution.tree.TreeUtils.isMonophyletic(
                (dr.evolution.tree.Tree) plain.objects.get("startingTree").getNativeObject(), pair));

        Run run = run(xml);
        assertTrue(dr.evolution.tree.TreeUtils.isMonophyletic(
                (dr.evolution.tree.Tree) run.objects.get("treeModel").getNativeObject(), pair));
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);
    }

    // estimated nucleotide frequencies are optimized on the simplex: stick-breaking transform, analytic gradient
    // (SubstitutionModelGradient behind branchSubstitutionParameterGradient), for HKY, TN93 and GTR, MAP and ML
    public void testFrequenciesOnTheSimplex() throws Exception {
        for (dr.evomodel.substmodel.nucleotide.NucModelType type : new dr.evomodel.substmodel.nucleotide.NucModelType[]{
                dr.evomodel.substmodel.nucleotide.NucModelType.HKY, dr.evomodel.substmodel.nucleotide.NucModelType.TN93,
                dr.evomodel.substmodel.nucleotide.NucModelType.GTR}) {
            for (MaximizationOptions.Objective objective : new MaximizationOptions.Objective[]{
                    MaximizationOptions.Objective.MAP, MaximizationOptions.Objective.ML}) {
                TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
                if (options == null) {
                    return;
                }
                for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                    model.setNucSubstitutionModel(type);
                }
                options.maximizationOptions.objective = objective;
                String xml = TreeLimeHarness.generate(quiet(options));
                String label = type + " " + objective;
                assertTrue(label, between(xml, "<multivariateCompoundTransform>", "</multivariateCompoundTransform>")
                        .contains("<UnitSimplexTransform dimension=\"4\"/>"));
                assertTrue(label, xml.contains("<branchSubstitutionParameterGradient id=\"gradient.frequencies\" homogeneous=\"true\">"));

                Run run = run(xml);
                assertEquals(label, MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
                assertTrue(label, run.end > run.start);
                Parameter frequencies = run.parameter("frequencies");
                double sum = 0;
                boolean moved = false;
                for (int i = 0; i < frequencies.getDimension(); i++) {
                    sum += frequencies.getParameterValue(i);
                    moved |= Math.abs(frequencies.getParameterValue(i) - 0.25) > 1E-3;
                }
                assertEquals(label, 1.0, sum, 1E-12);
                assertTrue(label + ": the frequencies moved", moved);
                assertGradientChecked(run);
            }
        }

        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        dr.app.treelime.options.Parameter frequencies = null;
        for (dr.app.treelime.options.Parameter parameter : options.selectParameters()) {
            if (parameter.getBaseName().equals("frequencies")) {
                frequencies = parameter;
            }
        }
        assertNotNull(frequencies);
        assertEquals(MaximizationOptions.Treatment.FREQUENCIES, options.maximizationOptions.classify(frequencies, false));
        // the all-numeric mode cannot keep them on the simplex
        options.maximizationOptions.setGradientMode(MaximizationOptions.GradientMode.NUMERIC);
        assertEquals(MaximizationOptions.Treatment.SIMPLEX_NUMERIC, options.maximizationOptions.classify(frequencies, false));
        options.maximizationOptions.setGradientMode(MaximizationOptions.GradientMode.ANALYTIC_WHERE_AVAILABLE);
        // a prior that is not flat on the simplex has no gradient term here
        frequencies.priorType = dr.app.treelime.types.PriorType.LOGNORMAL_PRIOR;
        assertEquals(MaximizationOptions.Treatment.SIMPLEX, options.maximizationOptions.classify(frequencies, false));
    }

    // a tree's node heights, the root included, are one row named nodeHeights; deselecting it holds them all fixed,
    // and files that deselected the root height still do
    public void testNodeHeightsRow() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        dr.app.treelime.options.Parameter rootHeight = null;
        for (dr.app.treelime.options.Parameter parameter : options.selectParameters()) {
            if (parameter.getName().equals("treeModel.rootHeight")) {
                rootHeight = parameter;
            }
        }
        assertNotNull(rootHeight);
        MaximizationOptions maximization = options.maximizationOptions;
        assertEquals("treeModel.nodeHeights", MaximizationOptions.displayName(rootHeight));
        assertEquals(MaximizationOptions.Treatment.NODE_HEIGHTS, maximization.classify(rootHeight, false));
        assertTrue(TreeLimeHarness.generate(quiet(options)).contains("<nodeHeightProxyParameter"));

        maximization.setExcluded(rootHeight, true);
        assertTrue(maximization.excluded.contains("treeModel.nodeHeights"));
        assertEquals(MaximizationOptions.Treatment.EXCLUDED, maximization.classify(rootHeight, false));
        assertFalse(TreeLimeHarness.generate(options).contains("<nodeHeightProxyParameter"));

        maximization.setExcluded(rootHeight, false);
        assertTrue(maximization.excluded.isEmpty());
        maximization.excluded.add("treeModel.rootHeight"); // as saved before the row was renamed
        assertEquals(MaximizationOptions.Treatment.EXCLUDED, maximization.classify(rootHeight, false));
        maximization.setExcluded(rootHeight, false);
        assertTrue(maximization.excluded.isEmpty());
    }

    // the fixed and random local clocks are not offered; a file that still has one is refused at Generate
    public void testLocalClocksAreRefused() throws Exception {
        for (ClockType type : new ClockType[]{ClockType.FIXED_LOCAL_CLOCK, ClockType.RANDOM_LOCAL_CLOCK}) {
            TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            assertFalse(type.isMaximizable());
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(type);
            }
            try {
                new dr.app.treelime.generator.BeastGenerator(options, TreeLimeHarness.COMPONENTS).checkOptions();
                fail("expected the " + type + " to be refused");
            } catch (dr.app.treelime.generator.Generator.GeneratorException e) {
                assertTrue(e.getMessage(), e.getMessage().contains("shrinkage local clock"));
                assertEquals(dr.app.treelime.TreeLimeFrame.CLOCK_MODELS, e.getSwitchToPanel());
            }
        }
        assertTrue(ClockType.SHRINKAGE_LOCAL_CLOCK.isMaximizable());
    }

    // Generate refuses a run time BEAST could not parse, and points at the Maximization tab
    public void testInvalidRunTimeIsRefused() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        options.maximizationOptions.maxRunTime = "5 lightyears";
        dr.app.treelime.generator.BeastGenerator generator =
                new dr.app.treelime.generator.BeastGenerator(options, TreeLimeHarness.COMPONENTS);
        try {
            generator.checkOptions();
            fail("expected the run time to be refused");
        } catch (dr.app.treelime.generator.Generator.GeneratorException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("lightyears"));
            assertEquals(dr.app.treelime.TreeLimeFrame.MAXIMIZATION, e.getSwitchToPanel());
        }
        options.maximizationOptions.maxRunTime = "30min";
        generator.checkOptions();
        assertTrue(TreeLimeHarness.generate(options).contains("maxRunTime=\"30min\""));
    }

    public void testMaximumLikelihood() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        options.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<objective>", "</objective>").contains("<likelihood idref=\"likelihood\"/>"));
        // the likelihood does not depend on the population size, so ML leaves it out; kappa has no prior term
        assertFalse(between(xml, "<compoundParameter id=\"optimizedParameters\">", "</compoundParameter>")
                .contains("constant.popSize"));
        assertFalse(xml.contains("<prior idref=\"prior\"/>\n\t\t\t<parameter idref=\"kappa\"/>"));

        Run run = run(xml);
        assertEquals("likelihood", run.maximizer.getLikelihood().getId());
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);
    }

    // a parameter unticked in the Maximization tab keeps its initial value
    public void testExcludedParameterStaysFixed() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        options.maximizationOptions.excluded.add("kappa");
        String xml = TreeLimeHarness.generate(quiet(options));
        assertFalse(between(xml, "<compoundParameter id=\"optimizedParameters\">", "</compoundParameter>")
                .contains("kappa"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertEquals(run.startValue("kappa"), run.parameter("kappa").getParameterValue(0), 0.0);
        assertTrue(run.parameter("constant.popSize").getParameterValue(0) != run.startValue("constant.popSize"));
    }

    public void testNothingToOptimize() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        options.maximizationOptions.excluded.add("kappa");
        options.maximizationOptions.excluded.add("constant.popSize");
        options.maximizationOptions.excluded.add("treeModel.rootHeight");
        options.maximizationOptions.excluded.add("frequencies");
        try {
            TreeLimeHarness.generate(quiet(options));
            fail("expected a generator exception");
        } catch (Exception e) {
            assertTrue(e.getMessage(), e.getMessage().contains("No parameter is selected"));
        }
    }

    // a vector parameter (the HMC clock's branch rates) whose dimension only BEAST knows: the transform finds its
    // segment from the compound parameter
    public void testVectorParameter() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            clock.setClockType(ClockType.HMC_CLOCK);
        }
        quiet(options).maximizationOptions.nIterations = 3;
        String xml = TreeLimeHarness.generate(options);
        assertTrue(between(xml, "<compoundParameter id=\"optimizedParameters\">", "</compoundParameter>")
                .contains("branchRatesHMC.rates"));

        assertTrue(between(xml, "<jointGradient id=\"gradient.branchRatesHMC.rates\">", "</jointGradient>")
                .contains("<branchRateGradient traitName=\"Sequence\">"));

        Run run = run(xml);
        Parameter rates = run.parameter("branchRatesHMC.rates");
        assertTrue(rates.getDimension() > 1);
        assertEquals(MaximizerWrtParameter.Reason.MAX_ITERATIONS, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);
    }

    // with tip dates BEAUti estimates the strict clock rate: the analytic branch-rate gradient
    public void testStrictClockRateWithTipDates() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
        if (options == null) {
            return;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        // with the node heights, rate and time trade off: L-BFGS needs more iterations
        quiet(options).maximizationOptions.nIterations = 2000;
        String xml = TreeLimeHarness.generate(options);
        assertTrue(between(xml, "<jointGradient id=\"gradient.clock.rate\">", "</jointGradient>")
                .contains("<branchRateGradient traitName=\"Sequence\">"));
        // the clock rate's CTMC reference prior depends on the tree length, so on the node heights too
        assertTrue(between(xml, "<jointGradient id=\"gradient.nodeHeights\">", "</jointGradient>")
                .contains("<ctmcScalePriorNodeHeightGradient>"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertTrue(run.parameter("clock.rate").getParameterValue(0) != 1.0);
        assertGradientChecked(run);
    }

    // the hmc_nodeHeight set-up: the root and internal heights as ratios on the real line, with the tree
    // likelihood's and the coalescent's analytic node-height gradients
    public void testNodeHeights() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(xml.contains("<nodeHeightProxyParameter id=\"nodeHeights\" rootNode=\"true\">"));
        assertTrue(xml.contains("<nodeHeightTransform id=\"nodeHeightTransform\" withRoot=\"true\">"));
        assertTrue(xml.contains("<ratios realLine=\"true\">"));
        String block = between(xml, "<jointGradient id=\"gradient.nodeHeights\">", "</jointGradient>");
        assertTrue(block.contains("<coalescentGradient>"));
        assertTrue(block.contains("<nodeHeightGradient traitName=\"Sequence\">"));
        String transform = between(xml, "<multivariateCompoundTransform>", "</multivariateCompoundTransform>");
        assertTrue(transform.indexOf("<compoundTransform>") < transform.indexOf("<nodeHeightTransform idref"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertTrue(run.parameter("treeModel.rootHeight").getParameterValue(0) != run.startValue("treeModel.rootHeight"));
        assertGradientChecked(run);

        // the prior side by finite differences, when a prior on the root height also depends on the heights
        options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        for (dr.app.treelime.options.Parameter parameter : options.selectParameters()) {
            if (parameter.getName().equals("treeModel.rootHeight")) {
                parameter.priorType = dr.app.treelime.types.PriorType.UNIFORM_PRIOR;
                parameter.uniformLower = 0.0;
                parameter.uniformUpper = 1E10;
            }
        }
        // BEAUti rescales the starting tree to the initial root height, which it computes as NaN without the GUI
        xml = TreeLimeHarness.generate(quiet(options)).replace("<rescaledTree id=\"startingTree\" height=\"NaN\">",
                "<rescaledTree id=\"startingTree\" height=\"30.0\">");
        block = between(xml, "<jointGradient id=\"gradient.nodeHeights\">", "</jointGradient>");
        assertTrue(block.contains("<numericalGradient>"));
        assertFalse(block.contains("<coalescentGradient>"));
        run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertGradientChecked(run);
    }

    // 20 states: the site-model gradient needs BEAGLE's transposed pre-order
    public void testAminoAcidsGamma() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/YFV.aa.fasta");
        if (options == null) {
            return;
        }
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setGammaHetero(true);
        }
        // 70 node heights: L-BFGS needs more iterations
        quiet(options).maximizationOptions.nIterations = 2000;
        String xml = TreeLimeHarness.generate(options);
        assertTrue(between(xml, "<jointGradient id=\"gradient.alpha\">", "</jointGradient>")
                .contains("<siteRateModelGradient>"));
        assertGammaShapeHeldAtBound(run(xml));
    }

    // with tip dates and free node heights the gamma shape of these data also runs to its boundary
    public void testDatedGammaShapeAtBoundary() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
        if (options == null) {
            return;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setGammaHetero(true);
        }
        quiet(options).maximizationOptions.nIterations = 2000;
        assertGammaShapeHeldAtBound(run(TreeLimeHarness.generate(options)));
    }

    // with free node heights the shape of these data runs to 0, where the discrete gamma tends to a model of invariant
    // sites and L-BFGS would chase a plateau until its line search fails: it is held at its default bound instead,
    // where the density still increases towards 0
    private static void assertGammaShapeHeldAtBound(Run run) {
        System.err.println(run.maximizer.getReport());
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertEquals(MaximizationOptions.GAMMA_SHAPE_LOWER_BOUND, run.parameter("alpha").getParameterValue(0), 0.0);
        assertEquals(1, run.maximizer.getHeld().size());
        assertTrue(run.maximizer.getHeld().get(0).toString(), run.maximizer.getHeld().get(0).gradient < 0);
        assertTrue(run.parameter("treeModel.rootHeight").getParameterValue(0) != run.startValue("treeModel.rootHeight"));
        // zero-length branches at the optimum put nodes next to tip dates, where the coalescent has kinks and central
        // differences along a random direction straddle them: check each coordinate at the optimum instead, against
        // whichever one-sided difference applies
        assertGradientMatchesAtOptimum(run);
    }

    private static void assertGradientMatchesAtOptimum(Run run) {
        com.github.lbfgs4j.liblbfgs.Function function = run.maximizer.getFunction();
        double[] y = run.maximizer.getMinimumPoint(false).clone();
        double[] gradient = function.gradientAt(y.clone());
        double value = function.valueAt(y.clone());
        final double h = 1E-6;
        for (int i = 0; i < y.length; i++) {
            if (gradient[i] == 0.0) {
                continue; // a held coordinate
            }
            double[] plus = y.clone();
            plus[i] += h;
            double[] minus = y.clone();
            minus[i] -= h;
            double up = function.valueAt(plus);
            double down = function.valueAt(minus);
            double tolerance = 1E-3 * Math.max(1.0, Math.abs(gradient[i]));
            double central = (up - down) / (2 * h);
            double forward = (up - value) / h;
            double backward = (value - down) / h;
            assertTrue("coordinate " + i + ": analytic " + gradient[i] + ", central " + central + ", forward " +
                            forward + ", backward " + backward,
                    Math.abs(gradient[i] - central) < tolerance || Math.abs(gradient[i] - forward) < tolerance ||
                            Math.abs(gradient[i] - backward) < tolerance);
        }
    }

    // each tree prior's analytic node-height gradient, joined with the tree likelihood's
    public void testNodeHeightTreePriors() throws Exception {
        String[][] cases = {
                {"EXPONENTIAL", "<coalescentGradient>"},
                {"SKYGRID", "<gmrfSkyrideGradient wrtParameter=\"nodeHeight\">"},
                {"GMRF_SKYRIDE", "<gmrfSkyrideGradient wrtParameter=\"nodeHeight\">"},
                {"YULE", "<speciationLikelihoodGradient wrtParameter=\"nodeHeight\">"},
                {"BIRTH_DEATH", "<speciationLikelihoodGradient wrtParameter=\"nodeHeight\">"},
        };
        for (String[] c : cases) {
            TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            for (dr.app.treelime.options.PartitionTreePrior prior : options.getPartitionTreePriors()) {
                prior.setNodeHeightPrior(dr.app.treelime.types.TreePriorType.valueOf(c[0]));
            }
            String xml = TreeLimeHarness.generate(quiet(options));
            String block = between(xml, "<jointGradient id=\"gradient.nodeHeights\">", "</jointGradient>");
            assertTrue(c[0] + " has " + c[1], block.contains(c[1]));
            assertTrue(block.contains("<nodeHeightGradient traitName=\"Sequence\">"));

            Run run = run(xml);
            System.err.println(c[0] + ": " + run.maximizer.getReason() + " after " + run.maximizer.getIterations() +
                    " iterations");
            assertTrue(c[0], run.end > run.start);
            assertGradientChecked(run);
        }
    }

    // finite differences are unreliable where the objective barely changes (a line-search trial that squeezes a node
    // onto its child, for example), so a rare large error is allowed; a missing or double-counted term fails
    // nearly every check
    private static void assertGradientChecked(Run run) {
        java.util.List<Double> errors = run.maximizer.getGradientCheckErrors();
        assertFalse(errors.isEmpty());
        int large = 0;
        for (double error : errors) {
            if (!(error < 1E-4)) {
                large++;
            }
        }
        assertTrue(large + " of " + errors.size() + " gradient checks have a relative error above 1E-4: " + errors,
                large <= Math.max(1, errors.size() / 20));
    }

    // the topology and node heights are fixed during maximization: start from UPGMA, not a random tree, so that
    // runs maximize the same model
    private static TreeLimeOptions quiet(TreeLimeOptions options) {
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            tree.setStartingTreeType(StartingTreeType.UPGMA);
        }
        options.maximizationOptions.nIterations = 200;
        options.maximizationOptions.printToScreen = false;
        return options;
    }

    private static String between(String xml, String from, String to) {
        int start = xml.indexOf(from);
        assertTrue("the XML has " + from, start >= 0);
        return xml.substring(start, xml.indexOf(to, start));
    }

    private static class Run {
        MaximizerWrtParameter maximizer;
        Map<String, XMLObject> objects;
        Map<String, Double> startValues = new java.util.HashMap<String, Double>();
        double start;
        double end;

        Parameter parameter(String id) {
            return (Parameter) objects.get(id).getNativeObject();
        }

        double startValue(String id) {
            return startValues.get(id);
        }
    }

    // parse without maximizing, so the start can be recorded, then maximize
    private static Run run(String xml) throws Exception {
        // check the gradient against finite differences of the objective at every evaluation
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\" checkGradientEvery=\"1\"");
        // the generated reports would write files and maximize nothing: drop them
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";

        Run run = new Run();
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        run.objects = parser.parse(new StringReader(xml), false);
        run.maximizer = (MaximizerWrtParameter) run.objects.get("maximizer").getNativeObject();
        for (String id : new String[]{"kappa", "constant.popSize", "treeModel.rootHeight"}) {
            if (run.objects.containsKey(id)) {
                run.startValues.put(id, run.parameter(id).getParameterValue(0));
            }
        }
        run.start = run.maximizer.getLikelihood().getLogLikelihood();
        run.maximizer.maximize();
        run.end = run.maximizer.getLikelihood().getLogLikelihood();
        return run;
    }
}
