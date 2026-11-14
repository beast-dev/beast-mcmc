/*
 * HierarchicalModelPruningTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.app.beauti;

import dr.app.beauti.components.hpm.HierarchicalModelComponentOptions;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.generator.Generator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.PriorType;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

/**
 * A hierarchical model holds its parameters directly, and changing the analysis does not destroy the
 * parameters it replaces, it only stops them being used. A model built beforehand could therefore be left
 * over parameters that are no longer in the analysis, and was written out regardless, giving a file that
 * referred to parameters it did not declare and that BEAST would not read.
 *
 * See <a href="https://github.com/beast-dev/beast-mcmc/issues/1062">issue #1062</a>.
 */
public class HierarchicalModelPruningTest {

    private static final String CASE = "apes_hky_strict_constant";

    private HierarchicalModelComponentOptions hpmOptions(BeautiOptions options) {
        return (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);
    }

    private Parameter active(BeautiOptions options, String name) {
        for (Parameter parameter : options.selectParameters()) {
            if (parameter.getName().equals(name)) {
                return parameter;
            }
        }
        return null;
    }

    /**
     * Partitions into three codon positions and puts the three kappas into one hierarchical model, as a
     * user would from the Priors panel.
     */
    private BeautiOptions withAModelOverTheCodonKappas(String name) throws Exception {
        BeautiOptions options = BeautiHeadlessHarness.getCase(CASE).build();
        BeautiHeadlessHarness.substitutionModel(options).setCodonHeteroPattern("123");

        List<Parameter> parameters = new ArrayList<Parameter>();
        for (String kappa : new String[]{"CP1.kappa", "CP2.kappa", "CP3.kappa"}) {
            Parameter parameter = active(options, kappa);
            assertNotNull("the test needs " + kappa + " to be in the analysis", parameter);
            parameter.priorType = PriorType.NORMAL_HPM_PRIOR;
            parameter.setPriorEdited(true);
            parameters.add(parameter);
        }
        hpmOptions(options).addHPM(name, parameters, PriorType.NORMAL_HPM_PRIOR);
        return options;
    }

    /** @return the message if generation was refused, or null if it succeeded */
    private String checkOptions(BeautiOptions options) {
        try {
            new BeastGenerator(options, BeautiHeadlessHarness.COMPONENTS).checkOptions();
            return null;
        } catch (Generator.GeneratorException e) {
            return e.getMessage();
        }
    }

    @Test
    public void aModelOverTheCurrentParametersIsLeftAlone() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");

        assertNull("nothing to prune", hpmOptions(options).pruneToActiveParameters(options.selectParameters()));
        assertNull("generation should not be refused", checkOptions(options));
        assertEquals(1, hpmOptions(options).getHPMList().size());
    }

    @Test
    public void aModelLosingAllButOneParameterIsRemoved() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");

        // 1+2 and 3 replaces CP1.kappa and CP2.kappa with CP1+2.kappa, leaving only CP3.kappa in the model
        BeautiHeadlessHarness.substitutionModel(options).setCodonHeteroPattern("112");

        String message = hpmOptions(options).pruneToActiveParameters(options.selectParameters());
        assertNotNull("the stale model should have been reported", message);
        assertTrue(message, message.contains("codonKappas"));
        assertTrue("the model should have been removed", hpmOptions(options).getHPMList().isEmpty());
    }

    @Test
    public void aParameterLeftBehindGetsAnOrdinaryPrior() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");
        BeautiHeadlessHarness.substitutionModel(options).setCodonHeteroPattern("112");

        String message = hpmOptions(options).pruneToActiveParameters(options.selectParameters());

        Parameter remaining = active(options, "CP3.kappa");
        assertNotNull("CP3.kappa is still in the analysis", remaining);
        assertFalse("it should no longer carry a hierarchical prior",
                remaining.priorType == PriorType.NORMAL_HPM_PRIOR
                        || remaining.priorType == PriorType.LOGNORMAL_HPM_PRIOR
                        || remaining.priorType == PriorType.EXPONENTIAL_HPM_PRIOR);
        assertNull("and should no longer be linked", remaining.linkedName);
        assertFalse(remaining.isLinked);
        assertTrue("the user should be told which parameter now has a default prior",
                message.contains("CP3.kappa"));
    }

    @Test
    public void aModelKeepingTwoParametersSurvivesWithoutTheRest() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");

        // only CP2.kappa leaves the analysis, so two members are left and the model is still a model
        List<Parameter> active = new ArrayList<Parameter>(options.selectParameters());
        active.remove(active(options, "CP2.kappa"));

        String message = hpmOptions(options).pruneToActiveParameters(active);
        assertNotNull(message);
        assertEquals("the model should survive", 1, hpmOptions(options).getHPMList().size());

        List<String> names = new ArrayList<String>();
        for (Parameter parameter : hpmOptions(options).getHPMList().get(0).getArgumentParameterList()) {
            names.add(parameter.getName());
        }
        assertEquals(2, names.size());
        assertTrue(names.contains("CP1.kappa"));
        assertTrue(names.contains("CP3.kappa"));
    }

    @Test
    public void generationIsRefusedUntilTheStaleModelHasGone() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");
        BeautiHeadlessHarness.substitutionModel(options).setCodonHeteroPattern("112");

        String refused = checkOptions(options);
        assertNotNull("the first attempt should be refused", refused);
        assertTrue(refused, refused.contains("codonKappas"));

        assertNull("the second attempt should succeed", checkOptions(options));
    }

    @Test
    public void theGeneratedFileDoesNotReferToParametersItDoesNotDeclare() throws Exception {
        BeautiOptions options = withAModelOverTheCodonKappas("codonKappas");
        BeautiHeadlessHarness.substitutionModel(options).setCodonHeteroPattern("112");
        checkOptions(options); // prunes, and is refused

        String xml = BeautiHeadlessHarness.generate("hpmPruned", options);

        assertFalse("the removed model should not be written", xml.contains("codonKappas"));
        assertFalse("nor should its parameters be referred to", xml.contains("idref=\"CP1.kappa\""));
        assertFalse(xml.contains("idref=\"CP2.kappa\""));
        assertTrue("the parameter still in the analysis should be declared",
                xml.contains("id=\"CP1+2.kappa\"") || xml.contains("id=\"CP3.kappa\""));
    }
}
