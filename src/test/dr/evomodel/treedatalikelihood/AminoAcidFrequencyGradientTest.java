/*
 * AminoAcidFrequencyGradientTest.java
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
import dr.evolution.datatype.AminoAcids;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.aminoacid.WAG;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient.Backend;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Analytic gradient wrt the frequencies of an empirical amino-acid model (WAG+F, 20 states; plan 5, FQ11(b, d)).
 */
public class AminoAcidFrequencyGradientTest extends TestCase {

    private static final String[][] SEQUENCES = {
            {"A", "ARNDCQEGHILKMFPSTWYVARNDCQEGHI"},
            {"B", "ARNDCQEGHILKMFPSTWYVSRNDCQEGHL"},
            {"C", "ARNECQEGHILRMFPSTWYVAKNDCQDGHI"},
            {"D", "SRNDCQEGHVLKMFPATWYVARNECQEGYI"},
            {"E", "ARQDCQEGHILKMFPSTFYVARNDCQEGHV"},
    };
    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final int SITE_COUNT = 30;

    private Parameter frequencies;
    private Parameter clockRate;
    private TreeDataLikelihood likelihood;

    public AminoAcidFrequencyGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(AminoAcids.INSTANCE);
        for (String[] row : SEQUENCES) {
            alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
        }
        Tree newick = new NewickImporter(NEWICK).importTree(null);
        DefaultTreeModel tree = new DefaultTreeModel("treeModel", newick);

        // skewed away from WAG's empirical frequencies so the gradient is not near zero
        double[] pi = new double[AminoAcids.INSTANCE.getStateCount()];
        double sum = 0.0;
        for (int k = 0; k < pi.length; k++) {
            pi[k] = 1.0 + k;
            sum += pi[k];
        }
        for (int k = 0; k < pi.length; k++) {
            pi[k] /= sum;
        }
        frequencies = new Parameter.Default(pi);
        clockRate = new Parameter.Default(1.0);

        EmpiricalAminoAcidModel wag = new EmpiricalAminoAcidModel(WAG.INSTANCE,
                new FrequencyModel(AminoAcids.INSTANCE, frequencies));
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                new SitePatterns(alignment, null, 0, -1, 1, true), new HomogeneousBranchModel(wag),
                new GammaSiteRateModel("siteModel", 0.5, 4), true, false, PartialsRescalingScheme.ALWAYS, false,
                new PreOrderSettings(true, false, true, true));
        likelihood = new TreeDataLikelihood(delegate, tree, new StrictClockBranchRates(clockRate));
    }

    public void testEnginesAgree() {
        double[] adjoint = gradient(Backend.ADJOINT).getGradientLogDensity();
        double[] perDimension = gradient(Backend.PER_DIMENSION).getGradientLogDensity();
        assertEquals(frequencies.getDimension(), adjoint.length);
        for (int k = 0; k < adjoint.length; k++) {
            assertEquals("k=" + k, adjoint[k], perDimension[k], 1E-8 * Math.max(1.0, Math.abs(adjoint[k])));
        }
    }

    public void testFrequencyWeightedSumMatchesClockDerivative() {
        // pi.g = N - r d log L / d r for any normalized Q
        double[] g = gradient(Backend.ADJOINT).getGradientLogDensity();
        double piG = 0.0;
        for (int k = 0; k < g.length; k++) {
            piG += frequencies.getParameterValue(k) * g[k];
        }

        final double h = 1E-6;
        clockRate.setParameterValue(0, 1.0 + h);
        double plus = likelihood.getLogLikelihood();
        clockRate.setParameterValue(0, 1.0 - h);
        double minus = likelihood.getLogLikelihood();
        clockRate.setParameterValue(0, 1.0);

        assertEquals(SITE_COUNT - (plus - minus) / (2 * h), piG, 1E-6 * Math.abs(piG));
    }

    public void testTangentCheckPasses() {
        String report = gradient(Backend.AUTO).getReport(); // throws if the tangent check fails
        assertTrue(report, report.contains("engine: adjoint"));
    }

    private HomogeneousSubstitutionModelGradient gradient(Backend backend) {
        return new HomogeneousSubstitutionModelGradient(likelihood, frequencies, backend, 1E-5);
    }
}
