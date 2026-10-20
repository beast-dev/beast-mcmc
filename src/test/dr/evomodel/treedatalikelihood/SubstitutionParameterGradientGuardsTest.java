/*
 * SubstitutionParameterGradientGuardsTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.DifferentialMassProvider.Mode;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionParameterGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.continuous.hmc.BranchSubstitutionParameterGradientParser;
import dr.inference.model.Parameter;
import dr.xml.XMLParseException;
import test.dr.inference.trace.TraceCorrelationAssert;

/**
 * Refusals on the per-dimension substitution-parameter gradient path (plan 5, FQ1 b and d).
 */
public class SubstitutionParameterGradientGuardsTest extends TraceCorrelationAssert {

    public SubstitutionParameterGradientGuardsTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testRefusesLikelihoodWithoutPreOrder() {
        Fixture f = new Fixture(PreOrderSettings.getDefault());
        assertPreOrderRefused(f);
    }

    public void testRefusesLikelihoodWithoutBranchInfinitesimalDerivative() {
        Fixture f = new Fixture(new PreOrderSettings(true, true, false, false));
        assertPreOrderRefused(f);
    }

    public void testAcceptsBranchInfinitesimalPreOrder() throws XMLParseException {
        Fixture f = new Fixture(branchInfinitesimal());
        assertSame(f.delegate, BranchSubstitutionParameterGradientParser.getPreOrderBeagleDelegate(f.likelihood));
    }

    public void testKappaGradientBuilds() {
        Fixture f = new Fixture(branchInfinitesimal());
        HomogeneousSubstitutionParameterGradient gradient = new HomogeneousSubstitutionParameterGradient(
                "kappa", f.likelihood, f.kappa, f.delegate, 0, Mode.EXACT, 1E-3);
        assertEquals(1, gradient.getGradientLogDensity().length);
    }

    public void testRefusesDimensionEqualToParameterDimension() {
        Fixture f = new Fixture(branchInfinitesimal());
        try {
            new HomogeneousSubstitutionParameterGradient("kappa", f.likelihood, f.kappa, f.delegate,
                    f.kappa.getDimension(), Mode.EXACT, 1E-3);
            fail("dim == parameter dimension must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    public void testRefusesDuplicateTraitName() {
        Fixture f = new Fixture(branchInfinitesimal());
        new HomogeneousSubstitutionParameterGradient("kappa", f.likelihood, f.kappa, f.delegate, 0, Mode.EXACT, 1E-3);
        try {
            new HomogeneousSubstitutionParameterGradient("kappa", f.likelihood, f.kappa, f.delegate, 0,
                    Mode.FIRST_ORDER, null);
            fail("a second gradient with the same traitName on one likelihood must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    public void testRefusesFrequencyParameter() {
        Fixture f = new Fixture(branchInfinitesimal());
        try {
            new HomogeneousSubstitutionParameterGradient("frequencies", f.likelihood, f.frequencies, f.delegate,
                    0, Mode.EXACT, 1E-3);
            fail("per-dimension frequency gradients must be refused (no root term, no simplex)");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    private static PreOrderSettings branchInfinitesimal() {
        return new PreOrderSettings(true, false, true, false);
    }

    private static void assertPreOrderRefused(Fixture f) {
        try {
            BranchSubstitutionParameterGradientParser.getPreOrderBeagleDelegate(f.likelihood);
            fail("a likelihood without usePreOrder and branchInfinitesimalDerivative must be refused");
        } catch (XMLParseException expected) {
            assertTrue(expected.getMessage().contains("branchInfinitesimalDerivative"));
        }
    }

    private class Fixture {

        final Parameter kappa = new Parameter.Default(2.0);
        final Parameter frequencies = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4});
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;

        Fixture(PreOrderSettings settings) {
            DefaultTreeModel tree = new DefaultTreeModel(ladderTree());
            HKY hky = new HKY(kappa, new FrequencyModel(Nucleotides.INSTANCE, frequencies));
            delegate = new BeagleDataLikelihoodDelegate(tree, new SitePatterns(alignment, null, 0, -1, 1, true),
                    new HomogeneousBranchModel(hky), new GammaSiteRateModel("siteModel"),
                    false, false, PartialsRescalingScheme.NONE, false, settings);
            likelihood = new TreeDataLikelihood(delegate, tree, new StrictClockBranchRates(new Parameter.Default(1.0)));
        }

        private SimpleTree ladderTree() {
            SimpleNode current = new SimpleNode();
            current.setTaxon(taxa[0]);
            current.setHeight(0.0);
            for (int i = 1; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                SimpleNode parent = new SimpleNode();
                parent.setHeight(0.05 * i);
                parent.addChild(current);
                parent.addChild(tip);
                current = parent;
            }
            return new SimpleTree(current);
        }
    }
}
