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
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * The XML that TreeLime's generator writes for the Maximization tab: it parses, the maximizer optimizes the chosen
 * density over the chosen parameters with a numeric gradient, and it stops by converging.
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
        assertTrue(xml.contains("<joint idref=\"joint\"/>"));
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
    }

    public void testMaximumLikelihood() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        options.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        String xml = TreeLimeHarness.generate(quiet(options));
        assertTrue(between(xml, "<densityWrtParameter>", "</densityWrtParameter>")
                .contains("<likelihood idref=\"likelihood\"/>"));

        Run run = run(xml);
        assertEquals("likelihood", run.maximizer.getLikelihood().getId());
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, run.maximizer.getReason());
        assertTrue(run.end > run.start);
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

        Run run = run(xml);
        Parameter rates = run.parameter("branchRatesHMC.rates");
        assertTrue(rates.getDimension() > 1);
        assertEquals(MaximizerWrtParameter.Reason.MAX_ITERATIONS, run.maximizer.getReason());
        assertTrue(run.end > run.start);
    }

    private static TreeLimeOptions quiet(TreeLimeOptions options) {
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
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
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
