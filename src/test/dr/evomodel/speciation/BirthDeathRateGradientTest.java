/*
 * BirthDeathRateGradientTest.java
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

package test.dr.evomodel.speciation;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.evomodel.speciation.SpeciationLikelihoodGradient.WrtParameter;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Collections;

/**
 * Gradients of the Gernhard (2008) birth-death tree prior wrt its parameters, as BEAUti writes them: the birth rate
 * r = lambda - mu (lambda for Yule), the relative death rate a = mu / lambda and the sampling proportion rho.
 */
public class BirthDeathRateGradientTest extends TestCase {

    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";

    public BirthDeathRateGradientTest(String name) {
        super(name);
    }

    public void testYule() throws Exception {
        Fixture f = new Fixture(2.0, Double.NaN, Double.NaN, false);
        check(f, WrtParameter.BIRTH_RATE, f.r);
    }

    public void testBirthDeath() throws Exception {
        for (double a : new double[]{0.1, 0.4, 0.9}) {
            Fixture f = new Fixture(1.5, a, Double.NaN, false);
            check(f, WrtParameter.BIRTH_RATE, f.r);
            check(f, WrtParameter.DEATH_RATE, f.a);
        }
    }

    public void testBirthDeathIncompleteSampling() throws Exception {
        for (double rho : new double[]{0.05, 0.3, 0.8}) {
            Fixture f = new Fixture(1.5, 0.4, rho, false);
            check(f, WrtParameter.BIRTH_RATE, f.r);
            check(f, WrtParameter.DEATH_RATE, f.a);
            check(f, WrtParameter.SAMPLING_PROBABILITY, f.rho);
        }
    }

    public void testConditionalOnRoot() throws Exception {
        Fixture f = new Fixture(1.5, 0.4, Double.NaN, true);
        check(f, WrtParameter.BIRTH_RATE, f.r);
        check(f, WrtParameter.DEATH_RATE, f.a);
    }

    public void testLargeRate() throws Exception {
        // e^{-r x} underflows to zero at every node (r x >= 1400 > 745); the gradient must stay finite and match
        Fixture f = new Fixture(20000.0, 0.4, 0.3, false);
        check(f, WrtParameter.BIRTH_RATE, f.r);
        check(f, WrtParameter.DEATH_RATE, f.a);
        check(f, WrtParameter.SAMPLING_PROBABILITY, f.rho);
    }

    public void testRefusals() throws Exception {
        Fixture yule = new Fixture(2.0, Double.NaN, Double.NaN, false);
        try {
            new SpeciationLikelihoodGradient(yule.likelihood, yule.tree, WrtParameter.DEATH_RATE);
            fail("a Yule model has no death rate");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("no relative death rate"));
        }

        SpeciationLikelihood excluding = new SpeciationLikelihood(yule.tree, yule.model,
                Collections.singleton(yule.tree.getNodeTaxon(yule.tree.getExternalNode(0))), "excluding");
        try {
            new SpeciationLikelihoodGradient(excluding, yule.tree, WrtParameter.BIRTH_RATE);
            fail("excluded taxa are not in the model's gradient");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("excluded taxa"));
        }
    }

    private static void check(Fixture f, WrtParameter wrt, Parameter parameter) {
        SpeciationLikelihoodGradient gradient = new SpeciationLikelihoodGradient(f.likelihood, f.tree, wrt);
        assertSame(parameter, gradient.getParameter());
        final double analytic = gradient.getGradientLogDensity()[0];
        assertTrue(Double.isFinite(analytic));

        final double value = parameter.getParameterValue(0);
        final double h = 1E-6 * value;
        parameter.setParameterValue(0, value + h);
        final double plus = f.likelihood.getLogLikelihood();
        parameter.setParameterValue(0, value - h);
        final double minus = f.likelihood.getLogLikelihood();
        parameter.setParameterValue(0, value);
        final double numeric = (plus - minus) / (2 * h);

        assertEquals(wrt + " at " + parameter.getId() + " = " + value, numeric, analytic,
                1E-6 * Math.max(1.0, Math.abs(numeric)));
    }

    private static class Fixture {

        final Parameter r;
        final Parameter a;
        final Parameter rho;
        final DefaultTreeModel tree;
        final BirthDeathGernhard08Model model;
        final SpeciationLikelihood likelihood;

        Fixture(double birth, double relativeDeath, double sampling, boolean conditionalOnRoot) throws Exception {
            Tree newick = new NewickImporter(NEWICK).importTree(null);
            tree = new DefaultTreeModel("treeModel", newick);
            r = new Parameter.Default("birthRate", birth);
            a = Double.isNaN(relativeDeath) ? null : new Parameter.Default("relativeDeathRate", relativeDeath);
            rho = Double.isNaN(sampling) ? null : new Parameter.Default("sampleProbability", sampling);
            model = new BirthDeathGernhard08Model("birthDeath", r, a, rho,
                    BirthDeathGernhard08Model.TreeType.UNSCALED, Units.Type.YEARS, conditionalOnRoot);
            likelihood = new SpeciationLikelihood(tree, model, "speciation");
        }
    }
}
