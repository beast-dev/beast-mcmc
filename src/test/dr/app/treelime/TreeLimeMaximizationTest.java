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

        // the all-numeric run reaches the same optimum
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

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);
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
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<jointGradient id=\"gradient.clock.rate\">", "</jointGradient>")
                .contains("<branchRateGradient traitName=\"Sequence\">"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertTrue(run.parameter("clock.rate").getParameterValue(0) != 1.0);
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
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<jointGradient id=\"gradient.alpha\">", "</jointGradient>")
                .contains("<siteRateModelGradient>"));

        Run run = run(xml);
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
        assertGradientChecked(run);
    }

    private static void assertGradientChecked(Run run) {
        String report = run.maximizer.getReport();
        String line = report.substring(report.indexOf("WorstGradientCheckError: "));
        double error = Double.parseDouble(line.substring(line.indexOf(' ') + 1, line.indexOf('\n')));
        assertTrue("worst relative gradient error " + error, error < 1E-4);
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
        for (String id : new String[]{"kappa", "constant.popSize"}) {
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
