/*
 * PartitionedGradientTest.java
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

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PartitionedGradients;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.SiteRateModelGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.hmc.CompoundParameterGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.ModelDependency;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Gradients across partitions, set up as BEAUti writes codon positions: three tree data likelihoods (sites 1, 2, 3
 * of the b1 alignment) sharing the tree, the strict clock, HKY (kappa, pi) and the gamma shape, with one relative rate
 * nu_i per partition (mu_i = nu_i W_i, W_i = 44 / n_i) collected in allNus, which sums to one.
 */
public class PartitionedGradientTest extends TestCase {

    private static final String[][] SEQUENCES = {
            {"A", "AAAAAACCCCGGGGGTTTTTTTAACAAAGTTACGARACCCGGTA"},
            {"B", "AAAAAACCCCGGGGGTTTTTTTAACGGCGTTATATAACCCGGCC"},
            {"C", "AAAAAACCCCGGGGGTTTTTTTGGTAAGGCCACGAA-CCCAATC"},
            {"D", "AAAAAACCCCGGGGGTTTTTTTGGTAATATTATATGGCCCGGCC"},
            {"E", "AAAAAACCCCGGGGGTTTTTTTGGTAAAATTGCGAGGTTTGGCA"},
    };
    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final int[] SITES = {15, 15, 14};

    private Parameter alpha;
    private Parameter kappa;
    private Parameter frequencies;
    private CompoundParameter allNus;
    private List<TreeDataLikelihood> partitions;
    private CompoundLikelihood likelihood;

    public PartitionedGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        for (String[] row : SEQUENCES) {
            alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
        }
        Tree newick = new NewickImporter(NEWICK).importTree(null);
        DefaultTreeModel tree = new DefaultTreeModel("treeModel", newick);
        StrictClockBranchRates clock = new StrictClockBranchRates(new Parameter.Default("clock.rate", 1.0));

        alpha = new Parameter.Default("alpha", 0.5);
        kappa = new Parameter.Default("kappa", 2.5);
        frequencies = new Parameter.Default("frequencies", new double[]{0.1, 0.2, 0.3, 0.4});
        HKY hky = new HKY(kappa, new FrequencyModel(Nucleotides.INSTANCE, frequencies));

        final double[] nu = {0.2, 0.3, 0.5};
        allNus = new CompoundParameter("allNus");
        partitions = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            Parameter nuI = new Parameter.Default("CP" + (i + 1) + ".nu", nu[i]);
            allNus.addParameter(nuI);
            GammaSiteRateModel siteModel = new GammaSiteRateModel("CP" + (i + 1) + ".siteModel", nuI,
                    44.0 / SITES[i], alpha, 4, null);
            SitePatterns patterns = new SitePatterns(alignment, null, i, -1, 3, true);
            BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree, patterns,
                    new HomogeneousBranchModel(hky), siteModel, true, false, PartialsRescalingScheme.NONE, false,
                    new PreOrderSettings(true, true, true, true));
            partitions.add(new TreeDataLikelihood(delegate, tree, clock));
        }
        List<Likelihood> members = new ArrayList<Likelihood>(partitions);
        likelihood = new CompoundLikelihood(members);
    }

    public void testAllNusBlockByBlock() {
        GradientWrtParameterProvider gradient = PartitionedGradients.create(partitions, allNus, SITE_RATE, 1.0, 1E-4);
        assertTrue(gradient instanceof CompoundParameterGradient);
        assertSame(allNus, gradient.getParameter());
        assertEquals(1.0, gradient.getSimplexSum());
        assertEquals(1.0, sum(allNus), 1E-15);

        // each block is nu_i of its own partition only
        assertGradient(gradient, allNus, 1E-6, 1E-7);
        gradient.getLikelihood(); // the compound of all three partitions
        assertEquals(3, gradient.getLikelihood().getLikelihoodSet().size());
        ((CompoundParameterGradient) gradient).getReport(); // tangent check on the simplex, throws on a mismatch
    }

    public void testLinkedShapeSummedOverPartitions() {
        GradientWrtParameterProvider gradient = PartitionedGradients.create(partitions, alpha, SITE_RATE,
                Double.NaN, 1E-4);
        assertTrue(gradient instanceof JointGradient);
        // the gamma quantiles make log L jitter at ~1e-9 in alpha, so differences check alpha to ~1e-4
        assertGradient(gradient, alpha, 1E-4, 1E-4);
    }

    public void testLinkedKappaAndFrequenciesSummedOverPartitions() {
        GradientWrtParameterProvider kappaGradient = PartitionedGradients.create(partitions, kappa, SUBSTITUTION,
                Double.NaN, 1E-5);
        assertTrue(kappaGradient instanceof JointGradient);
        assertGradient(kappaGradient, kappa, 1E-6, 1E-7);

        GradientWrtParameterProvider frequencyGradient = PartitionedGradients.create(partitions, frequencies,
                SUBSTITUTION, Double.NaN, 1E-5);
        assertEquals(1.0, frequencyGradient.getSimplexSum());
        ((JointGradient) frequencyGradient).getReport(); // tangent check, throws on a mismatch
    }

    public void testSinglePartitionIsUnwrapped() {
        List<TreeDataLikelihood> one = partitions.subList(0, 1);
        assertTrue(PartitionedGradients.create(one, allNus.getParameter(0), SITE_RATE, Double.NaN, 1E-4)
                instanceof SiteRateModelGradient);
    }

    public void testRefusals() {
        try {
            PartitionedGradients.create(partitions, new Parameter.Default("unused", 1.0), SITE_RATE, Double.NaN, 1E-4);
            fail("a parameter no partition uses must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not used by any partition"));
        }
        try {
            PartitionedGradients.create(partitions, alpha, SITE_RATE, 1.0, 1E-4);
            fail("simplexSum on a parameter differentiated whole must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("simplexSum"));
        }
        try {
            PartitionedGradients.treeDataLikelihoods(new CompoundLikelihood(new ArrayList<Likelihood>() {{
                add(partitions.get(0));
                add(new dr.inference.model.DummyLikelihood(null));
            }}));
            fail("a compound with a non-tree-data likelihood must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not a treeDataLikelihood"));
        }
    }

    private static final PartitionedGradients.Factory SITE_RATE = new PartitionedGradients.Factory() {
        @Override
        public boolean appliesTo(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
            return SiteRateModelGradient.isSiteRateParameter(treeDataLikelihood, parameter);
        }

        @Override
        public GradientWrtParameterProvider create(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
            return new SiteRateModelGradient(treeDataLikelihood, parameter, 1E-4);
        }
    };

    private static final PartitionedGradients.Factory SUBSTITUTION = new PartitionedGradients.Factory() {
        @Override
        public boolean appliesTo(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
            return !ModelDependency.findHolders(treeDataLikelihood, parameter).isEmpty();
        }

        @Override
        public GradientWrtParameterProvider create(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
            return new HomogeneousSubstitutionModelGradient(treeDataLikelihood, parameter,
                    HomogeneousSubstitutionModelGradient.Backend.ADJOINT, 1E-5);
        }
    };

    private void assertGradient(GradientWrtParameterProvider gradient, Parameter parameter, double relativeStep,
                                double tolerance) {
        double[] analytic = gradient.getGradientLogDensity();
        assertEquals(parameter.getDimension(), analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            final double value = parameter.getParameterValue(i);
            final double h = relativeStep * Math.max(1.0, value);
            parameter.setParameterValue(i, value + h);
            final double plus = likelihood.getLogLikelihood();
            parameter.setParameterValue(i, value - h);
            final double minus = likelihood.getLogLikelihood();
            parameter.setParameterValue(i, value);
            final double numeric = (plus - minus) / (2 * h);
            assertEquals(parameter.getId() + "[" + i + "]", numeric, analytic[i],
                    tolerance * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static double sum(Parameter parameter) {
        double sum = 0.0;
        for (int i = 0; i < parameter.getDimension(); i++) {
            sum += parameter.getParameterValue(i);
        }
        return sum;
    }
}
