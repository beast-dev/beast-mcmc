/*
 * HierarchicalRepresentationTest.java
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

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.InstanceDetails;
import beagle.ResourceDetails;
import dr.evolution.alignment.HierarchicalRateAggregatedPattern;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.AbstractHierarchicalRateProvider;
import dr.evomodel.substmodel.EqualRatesProvider;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.RelativeRatesProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BeagleFunctionality;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.substmodel.HierarchicalSubstitutionModelParser;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

/**
 * @author Xiang Ji
 */
public class HierarchicalRepresentationTest extends TraceCorrelationAssert {

    private static final int S = 4;

    private TreeModel treeModel;
    private final NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public HierarchicalRepresentationTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        format.setMaximumFractionDigits(5);
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();
    }

    private static boolean hierarchicalAvailable() {
        try {
            if (!BeagleFunctionality.IS_HIERARCHICAL_SUPPORTED()) {
                return false;
            }
            for (ResourceDetails resource : BeagleFactory.getResourceDetails()) {
                if ((resource.getFlags() & BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()) != 0) {
                    return true;
                }
            }
        } catch (Throwable ignored) {
        }
        return false;
    }

    private static boolean skip(String test) {
        if (hierarchicalAvailable()) {
            return false;
        }
        System.out.println("Skipping " + test + ": no BEAGLE resource provides the hierarchical representation");
        return true;
    }

    private static String clade(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return tree.getTaxonId(node.getNumber());
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getChildCount(node); i++) {
            sb.append(clade(tree, tree.getChild(node, i)));
        }
        char[] chars = sb.toString().toCharArray();
        Arrays.sort(chars);
        return new String(chars);
    }

    private static HierarchicalSubstitutionModel hierarchical(String newick, double[] pi, double lambdaRoot,
                                                              double lambdaAG, double lambdaCT) throws Exception {
        Tree imported = new NewickImporter(newick).importTree(null);
        TreeModel characterTree = new DefaultTreeModel("characterTree", imported);
        Parameter rates = new Parameter.Default(characterTree.getInternalNodeCount(), 0.0);
        for (int i = 0; i < characterTree.getInternalNodeCount(); i++) {
            String clade = clade(characterTree, characterTree.getNode(S + i));
            rates.setParameterValue(i, clade.equals("AG") ? lambdaAG : clade.equals("CT") ? lambdaCT : lambdaRoot);
        }
        return new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, characterTree),
                new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(pi)));
    }

    private TreeDataLikelihood likelihood(SubstitutionModel model, GammaSiteRateModel siteRateModel,
                                          boolean useHierarchicalRepresentation) {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                new HomogeneousBranchModel(model),
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.DEFAULT,
                false,
                new PreOrderSettings(false, false, false, true, false, false, false, useHierarchicalRepresentation));
        return new TreeDataLikelihood(delegate, treeModel, new DefaultBranchRateModel());
    }

    public void testTransitionMatricesFromBeagle() throws Exception {
        if (skip("testTransitionMatricesFromBeagle")) {
            return;
        }
        final double[] pi = {0.1, 0.2, 0.3, 0.4};
        final double[] categoryRates = {0.3, 1.0, 2.7};
        final int C = categoryRates.length;

        Beagle beagle = BeagleFactory.loadBeagleInstance(3, 5, 0, S, 10, 1, 1, C, 0, null,
                BeagleFlag.PRECISION_DOUBLE.getMask(),
                BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask() | BeagleFlag.PRECISION_DOUBLE.getMask());
        InstanceDetails details = beagle.getDetails();
        assertTrue("hierarchical instance", details != null &&
                (details.getFlags() & BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()) != 0);
        beagle.setStateFrequencies(0, pi);

        for (String newick : new String[]{"((A:1,G:1):1,(C:1,T:1):1);", "((T:1,C:1):1,(G:1,A:1):1);"}) {
            HierarchicalSubstitutionModel model = hierarchical(newick, pi, 1.3, 0.7, 2.1);
            double[] bundle = new double[model.getBundleLength(C)];
            double[] fromBeagle = new double[C * S * S];
            double[] fromJava = new double[S * S];
            for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
                model.fillBundle(t, categoryRates, bundle);
                beagle.setTransitionMatrix(0, bundle, 0.0);
                beagle.getTransitionMatrix(0, fromBeagle);
                for (int l = 0; l < C; l++) {
                    model.getTransitionProbabilities(t * categoryRates[l], fromJava);
                    for (int k = 0; k < S * S; k++) {
                        assertEquals(newick + " t=" + t + " category " + l + " entry " + k,
                                fromJava[k], fromBeagle[l * S * S + k], 1e-14);
                    }
                }
            }
        }
    }

    private void checkAgainstPublished(String name, double[] pi, GammaSiteRateModel hierarchicalSites,
                                       GammaSiteRateModel referenceSites, double expected) throws Exception {
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(pi));
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 8.0);
        HKY hky = new HKY(kappa, frequencies);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                RelativeRatesProvider.hky(kappa, frequencies));

        double hierarchical = likelihood(model, hierarchicalSites, true).getLogLikelihood();
        double eigen = likelihood(model, hierarchicalSites, false).getLogLikelihood();
        double reference = likelihood(hky, referenceSites, false).getLogLikelihood();

        assertEquals(name + " published", format.format(expected), format.format(hierarchical));
        assertEquals(name + " vs eigen", eigen, hierarchical, 1e-10 * Math.abs(eigen));
        assertEquals(name + " vs HKY", reference, hierarchical, 1e-10 * Math.abs(reference));
    }

    private static GammaSiteRateModel gamma(double alpha, int categories) {
        return new GammaSiteRateModel("gamma", null, 1.0, new Parameter.Default(alpha), categories, null);
    }

    public void testLikelihoodK80() throws Exception {
        if (skip("testLikelihoodK80")) {
            return;
        }
        double[] pi = {0.25, 0.25, 0.25, 0.25};
        checkAgainstPublished("K80", pi, new GammaSiteRateModel("gamma"), new GammaSiteRateModel("gamma"), -1868.89782);
    }

    public void testLikelihoodHKY85() throws Exception {
        if (skip("testLikelihoodHKY85")) {
            return;
        }
        double[] pi = alignment.getStateFrequencies();
        checkAgainstPublished("HKY85", pi, new GammaSiteRateModel("gamma"), new GammaSiteRateModel("gamma"), -1839.84514);
    }

    public void testLikelihoodHKY85G() throws Exception {
        if (skip("testLikelihoodHKY85G")) {
            return;
        }
        double[] pi = alignment.getStateFrequencies();
        checkAgainstPublished("HKY85G", pi, gamma(0.5, 4), gamma(0.5, 4), -1816.82611);
    }

    public void testTN93ThroughBeagle() throws Exception {
        if (skip("testTN93ThroughBeagle")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        Parameter kappa1 = new Parameter.Default(1.7), kappa2 = new Parameter.Default(3.1);
        TN93 tn93 = new TN93(kappa1, kappa2, frequencies);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                RelativeRatesProvider.tn93(kappa1, kappa2, frequencies));
        double hierarchical = likelihood(model, gamma(0.5, 4), true).getLogLikelihood();
        double reference = likelihood(tn93, gamma(0.5, 4), false).getLogLikelihood();
        assertEquals("TN93 vs dense", reference, hierarchical, 1e-10 * Math.abs(reference));
    }

    public void testKappaChangeAndRestoreThroughBeagle() throws Exception {
        if (skip("testKappaChangeAndRestoreThroughBeagle")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 8.0);
        HKY hky = new HKY(kappa, frequencies);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                RelativeRatesProvider.hky(kappa, frequencies));
        TreeDataLikelihood hierarchical = likelihood(model, new GammaSiteRateModel("gamma"), true);
        TreeDataLikelihood reference = likelihood(hky, new GammaSiteRateModel("gamma"), false);

        final double first = hierarchical.getLogLikelihood();
        assertEquals("kappa 8", reference.getLogLikelihood(), first, 1e-10 * Math.abs(first));
        hierarchical.storeModelState();
        reference.storeModelState();
        kappa.setParameterValue(0, 3.0);
        final double moved = hierarchical.getLogLikelihood();
        assertEquals("kappa 3", reference.getLogLikelihood(), moved, 1e-10 * Math.abs(moved));
        assertTrue("the likelihood moved", Math.abs(moved - first) > 1e-6);
        hierarchical.restoreModelState();
        reference.restoreModelState();
        assertEquals("restored", first, hierarchical.getLogLikelihood(), 0.0);
        kappa.setParameterValue(0, 5.0);
        assertEquals("after restore", reference.getLogLikelihood(), hierarchical.getLogLikelihood(),
                1e-10 * Math.abs(first));
    }

    public void testNegativeRatesRejected() throws Exception {
        if (skip("testNegativeRatesRejected")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 0.5);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                RelativeRatesProvider.hky(kappa, frequencies));
        assertFalse(model.hasNonNegativeRates());
        try {
            likelihood(model, new GammaSiteRateModel("gamma"), true);
            fail("kappa 0.5 must be refused at construction");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("clade AG has rate -")
                    && e.getMessage().contains("non-negative clade rates"));
        }
        assertTrue("the eigen path keeps kappa < 1",
                !Double.isInfinite(likelihood(model, new GammaSiteRateModel("gamma"), false).getLogLikelihood()));

        kappa.setParameterValue(0, 2.0);
        HKY hky = new HKY(kappa, frequencies);
        TreeDataLikelihood hierarchical = likelihood(model, new GammaSiteRateModel("gamma"), true);
        TreeDataLikelihood reference = likelihood(hky, new GammaSiteRateModel("gamma"), false);
        final double first = hierarchical.getLogLikelihood();
        assertEquals(reference.getLogLikelihood(), first, 1e-10 * Math.abs(first));

        hierarchical.storeModelState();
        kappa.setParameterValue(0, 0.5);
        assertEquals("kappa 0.5 rejected", Double.NEGATIVE_INFINITY, hierarchical.getLogLikelihood());
        hierarchical.restoreModelState();
        assertEquals("restored", first, hierarchical.getLogLikelihood(), 0.0);

        kappa.setParameterValue(0, 0.5);
        assertEquals(Double.NEGATIVE_INFINITY, hierarchical.getLogLikelihood());
        kappa.setParameterValue(0, 3.0);
        assertEquals("valid again", reference.getLogLikelihood(), hierarchical.getLogLikelihood(),
                1e-10 * Math.abs(first));
    }

    public void testEigenPathInsideTheGeneratorCone() throws Exception {
        if (skip("testEigenPathInsideTheGeneratorCone")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 0.5);
        HKY hky = new HKY(kappa, frequencies);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                RelativeRatesProvider.hky(kappa, frequencies));
        double eigen = likelihood(model, gamma(0.5, 4), false).getLogLikelihood();
        double reference = likelihood(hky, gamma(0.5, 4), false).getLogLikelihood();
        assertEquals("kappa 0.5, negative clade rates, eigen path", reference, eigen, 1e-10 * Math.abs(reference));
    }

    public void testNegativeRateInsideTheConeOnBothPaths() throws Exception {
        if (skip("testNegativeRateInsideTheConeOnBothPaths")) {
            return;
        }
        Parameter pi = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4});
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE, pi);
        TreeModel characterTree = new DefaultTreeModel("characterTree",
                new NewickImporter("((A:1,G:1):1,(C:1,T:1):1);").importTree(null));
        Parameter rates = new Parameter.Default(3, 0.0);
        for (int i = 0; i < 3; i++) {
            String clade = clade(characterTree, characterTree.getNode(S + i));
            rates.setParameterValue(i, clade.equals("AG") ? -0.35 : clade.equals("CT") ? 0.2 : 1.0);
        }
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, characterTree, frequencies, false));
        assertTrue(model.isValidGenerator());
        assertTrue(rates.isWithinBounds());
        TreeDataLikelihood eigen = likelihood(model, new GammaSiteRateModel("gamma"), false);
        final double first = eigen.getLogLikelihood();
        assertFalse(Double.isInfinite(first));
        try {
            likelihood(model, new GammaSiteRateModel("gamma"), true);
            fail("the plugin cannot take the negative AG rate");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("the generator is valid"));
        }
        pi.setParameterValue(1, 0.3);
        pi.setParameterValue(2, 0.2);
        assertFalse("pi_AG 0.3 puts sigma(AG) below 0", model.isValidGenerator());
    }

    public void testLikelihoodJC69EqualRates() throws Exception {
        if (skip("testLikelihoodJC69EqualRates")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25}));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 1.0), frequencies);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                new EqualRatesProvider(AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE),
                        frequencies, true));
        double hierarchical = likelihood(model, new GammaSiteRateModel("gamma"), true).getLogLikelihood();
        double reference = likelihood(hky, new GammaSiteRateModel("gamma"), false).getLogLikelihood();
        assertEquals("JC69 published", format.format(-1992.20564), format.format(hierarchical));
        assertEquals("JC69 vs HKY kappa 1", reference, hierarchical, 1e-10 * Math.abs(reference));
    }

    private void checkBundlesThroughBeagle(String what, HierarchicalSubstitutionModel model) {
        final int n = model.getDataType().getStateCount();
        final double[] categoryRates = {0.3, 1.0, 2.7};
        final int C = categoryRates.length;
        Beagle beagle = BeagleFactory.loadBeagleInstance(3, 5, 0, n, 10, 1, 1, C, 0, null,
                BeagleFlag.PRECISION_DOUBLE.getMask(),
                BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask() | BeagleFlag.PRECISION_DOUBLE.getMask());
        InstanceDetails details = beagle.getDetails();
        assertTrue(what + " hierarchical instance", details != null &&
                (details.getFlags() & BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()) != 0);
        beagle.setStateFrequencies(0, model.getFrequencyModel().getFrequencies());
        double[] bundle = new double[model.getBundleLength(C)];
        double[] fromBeagle = new double[C * n * n];
        double[] fromJava = new double[n * n];
        for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
            model.fillBundle(t, categoryRates, bundle);
            beagle.setTransitionMatrix(0, bundle, 0.0);
            beagle.getTransitionMatrix(0, fromBeagle);
            for (int l = 0; l < C; l++) {
                model.getTransitionProbabilities(t * categoryRates[l], fromJava);
                for (int k = 0; k < n * n; k++) {
                    assertEquals(what + " t=" + t + " category " + l + " entry " + k,
                            fromJava[k], fromBeagle[l * n * n + k], 1e-14);
                }
            }
        }
    }

    public void testEqualRatesBundlesThroughBeagle() throws Exception {
        if (skip("testEqualRatesBundlesThroughBeagle")) {
            return;
        }
        FrequencyModel binary = new FrequencyModel(dr.evolution.datatype.TwoStates.INSTANCE,
                new Parameter.Default(new double[]{0.3, 0.7}));
        checkBundlesThroughBeagle("two states", new HierarchicalSubstitutionModel("hierarchical",
                new EqualRatesProvider(AbstractHierarchicalRateProvider.defaultCharacterTree(binary.getDataType()),
                        binary, true)));

        dr.evolution.datatype.DataType seven = new dr.evolution.datatype.GeneralDataType(
                new String[]{"a", "b", "c", "d", "e", "f", "g"});
        FrequencyModel sevenFrequencies = new FrequencyModel(seven,
                new Parameter.Default(new double[]{0.05, 0.1, 0.12, 0.13, 0.15, 0.2, 0.25}));
        checkBundlesThroughBeagle("seven states", new HierarchicalSubstitutionModel("hierarchical",
                new EqualRatesProvider(AbstractHierarchicalRateProvider.defaultCharacterTree(seven), sevenFrequencies, true)));
    }

    public void testDefaultTopologyChangeThroughBeagle() throws Exception {
        if (skip("testDefaultTopologyChangeThroughBeagle")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        TreeModel characterTree = new DefaultTreeModel("characterTree",
                new NewickImporter("((A:1,G:1):1,(C:1,T:1):1);").importTree(null));
        double[] values = {0.6, 2.1, 0.9};
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(values.clone()),
                        characterTree, frequencies, false));
        TreeDataLikelihood hierarchical = likelihood(model, new GammaSiteRateModel("gamma"), true);
        TreeDataLikelihood eigen = likelihood(model, new GammaSiteRateModel("gamma"), false);
        final double first = hierarchical.getLogLikelihood();
        assertEquals(eigen.getLogLikelihood(), first, 1e-10 * Math.abs(first));

        hierarchical.storeModelState();
        eigen.storeModelState();
        dr.evolution.tree.NodeRef g = null, c = null;
        for (int i = 0; i < characterTree.getExternalNodeCount(); i++) {
            if (characterTree.getTaxonId(i).equals("G")) g = characterTree.getExternalNode(i);
            if (characterTree.getTaxonId(i).equals("C")) c = characterTree.getExternalNode(i);
        }
        dr.evolution.tree.NodeRef pg = characterTree.getParent(g), pc = characterTree.getParent(c);
        characterTree.beginTreeEdit();
        characterTree.removeChild(pg, g);
        characterTree.removeChild(pc, c);
        characterTree.addChild(pg, c);
        characterTree.addChild(pc, g);
        characterTree.endTreeEdit();

        HierarchicalSubstitutionModel fresh = new HierarchicalSubstitutionModel("fresh",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(values.clone()),
                        characterTree, frequencies, false));
        final double expected = likelihood(fresh, new GammaSiteRateModel("gamma"), false).getLogLikelihood();
        final double moved = hierarchical.getLogLikelihood();
        assertTrue("the topology change matters", Math.abs(moved - first) > 1e-3);
        assertEquals("hierarchical path after the change", expected, moved, 1e-10 * Math.abs(expected));
        assertEquals("eigen path after the change", expected, eigen.getLogLikelihood(), 1e-10 * Math.abs(expected));

        hierarchical.restoreModelState();
        eigen.restoreModelState();
        assertEquals("hierarchical path restored", first, hierarchical.getLogLikelihood(), 0.0);
        assertEquals("eigen path restored", first, eigen.getLogLikelihood(), 1e-10 * Math.abs(first));
    }

    private static String balancedNewick(String[] codes, int from, int to) {
        if (to - from == 1) {
            return codes[from] + ":1";
        }
        int mid = (from + to) / 2;
        return "(" + balancedNewick(codes, from, mid) + "," + balancedNewick(codes, mid, to) + "):1";
    }

    private static HierarchicalSubstitutionModel generalModel(String name, TreeModel characterTree, double[] rates,
                                                              FrequencyModel frequencies) {
        return new HierarchicalSubstitutionModel(name,
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(rates.clone()),
                        characterTree, frequencies, false));
    }

    private static void swapTips(TreeModel characterTree, java.util.Random random) {
        final int tips = characterTree.getExternalNodeCount();
        NodeRef a, b;
        do {
            a = characterTree.getExternalNode(random.nextInt(tips));
            b = characterTree.getExternalNode(random.nextInt(tips));
        } while (characterTree.getParent(a) == characterTree.getParent(b));
        NodeRef pa = characterTree.getParent(a), pb = characterTree.getParent(b);
        characterTree.beginTreeEdit();
        characterTree.removeChild(pa, a);
        characterTree.removeChild(pb, b);
        characterTree.addChild(pa, b);
        characterTree.addChild(pb, a);
        characterTree.endTreeEdit();
    }

    public void testCharacterTreeMovesWithSparseStatesThroughBeagle() throws Exception {
        if (skip("testCharacterTreeMovesWithSparseStatesThroughBeagle")) {
            return;
        }
        final int states = 16, columns = 150;
        String[] codes = new String[states];
        for (int i = 0; i < states; i++) {
            codes[i] = String.valueOf((char) ('A' + i));
        }
        dr.evolution.datatype.GeneralDataType general = new dr.evolution.datatype.GeneralDataType(codes);
        java.util.Random random = new java.util.Random(20261103);
        String[][] taxonSequence = new String[2][];
        taxonSequence[0] = PRIMATES_TAXON_SEQUENCE[0];
        taxonSequence[1] = new String[taxonSequence[0].length];
        StringBuilder[] rows = new StringBuilder[taxonSequence[0].length];
        for (int t = 0; t < rows.length; t++) {
            rows[t] = new StringBuilder();
        }
        for (int k = 0; k < columns; k++) {
            int[] pool = new int[1 + random.nextInt(3)];
            for (int j = 0; j < pool.length; j++) {
                pool[j] = random.nextInt(states);
            }
            for (StringBuilder row : rows) {
                row.append(codes[pool[random.nextInt(pool.length)]]);
            }
        }
        for (int t = 0; t < rows.length; t++) {
            taxonSequence[1][t] = rows[t].toString();
        }
        createAlignment(taxonSequence, general);
        treeModel = createPrimateTreeModel();
        treeModel.beginTreeEdit();
        for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {
            NodeRef node = treeModel.getInternalNode(i);
            treeModel.setNodeHeight(node, 20.0 * treeModel.getNodeHeight(node));
        }
        treeModel.endTreeEdit();

        double[] pi = new double[states];
        double total = 0.0;
        for (int i = 0; i < states; i++) {
            pi[i] = 0.5 + random.nextDouble();
            total += pi[i];
        }
        for (int i = 0; i < states; i++) {
            pi[i] /= total;
        }
        FrequencyModel frequencies = new FrequencyModel(general, new Parameter.Default(pi));
        TreeModel characterTree = new DefaultTreeModel("characterTree",
                new NewickImporter(balancedNewick(codes, 0, states) + ";").importTree(null));
        double[] rates = new double[states - 1];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = 0.3 + 1.7 * random.nextDouble();
        }
        Parameter rateParameter = new Parameter.Default(rates.clone());
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rateParameter, characterTree,
                        frequencies, false));

        for (int categories : new int[]{1, 4}) {
            GammaSiteRateModel sites = categories == 1 ? new GammaSiteRateModel("gamma") : gamma(0.5, categories);
            TreeDataLikelihood hierarchical = likelihood(model, sites, true);
            TreeDataLikelihood eigen = likelihood(model, sites, false);
            double current = hierarchical.getLogLikelihood();
            assertEquals("C=" + categories + " start", eigen.getLogLikelihood(), current, 1e-10 * Math.abs(current));
            int topologyMoves = 0, rejected = 0;
            for (int step = 0; step < 18; step++) {
                hierarchical.storeModelState();
                eigen.storeModelState();
                final int kind = step % 3;
                if (kind == 0) {
                    swapTips(characterTree, random);
                    topologyMoves++;
                } else if (kind == 1) {
                    int i = random.nextInt(rateParameter.getDimension());
                    rateParameter.setParameterValue(i, rateParameter.getParameterValue(i) * (0.7 + 0.6 * random.nextDouble()));
                } else {
                    NodeRef node = treeModel.getInternalNode(random.nextInt(treeModel.getInternalNodeCount()));
                    if (!treeModel.isRoot(node)) {
                        double low = Math.max(treeModel.getNodeHeight(treeModel.getChild(node, 0)),
                                treeModel.getNodeHeight(treeModel.getChild(node, 1)));
                        double high = treeModel.getNodeHeight(treeModel.getParent(node));
                        treeModel.setNodeHeight(node, low + (high - low) * random.nextDouble());
                    }
                }
                final String what = "C=" + categories + " step " + step;
                final double proposed = hierarchical.getLogLikelihood();
                final double fresh = likelihood(generalModel("fresh", characterTree, rateParameter.getParameterValues(),
                        frequencies), sites, false).getLogLikelihood();
                assertEquals(what + " hierarchical against a fresh eigen model", fresh, proposed, 1e-10 * Math.abs(fresh));
                assertEquals(what + " eigen path", fresh, eigen.getLogLikelihood(), 1e-10 * Math.abs(fresh));
                if ((step / 3) % 2 == 1) {
                    hierarchical.restoreModelState();
                    eigen.restoreModelState();
                    rejected++;
                    assertEquals(what + " restored", current, hierarchical.getLogLikelihood(), 0.0);
                } else {
                    hierarchical.acceptModelState();
                    eigen.acceptModelState();
                    current = proposed;
                }
            }
            assertEquals(6, topologyMoves);
            assertEquals(9, rejected);
            hierarchical.makeDirty();
            eigen.makeDirty();
            final double fresh = likelihood(generalModel("fresh", characterTree, rateParameter.getParameterValues(),
                    frequencies), sites, false).getLogLikelihood();
            assertEquals("C=" + categories + " full recompute at the end", fresh, hierarchical.getLogLikelihood(),
                    1e-10 * Math.abs(fresh));
            assertEquals("C=" + categories + " eigen full recompute at the end", fresh, eigen.getLogLikelihood(),
                    1e-10 * Math.abs(fresh));
        }
    }

    // Step 9 of the aggregation plan: the aggregate attribute reaches the plugin as a flag bit.

    private static final HierarchicalSubstitutionModel.Aggregation[] AGGREGATIONS = {
            HierarchicalSubstitutionModel.Aggregation.ON, HierarchicalSubstitutionModel.Aggregation.OFF,
            HierarchicalSubstitutionModel.Aggregation.AUTO};

    private static String randomNewick(java.util.List<String> taxa, java.util.Random random) {
        if (taxa.size() == 1) {
            return taxa.get(0) + ":" + (0.02 + 0.2 * random.nextDouble());
        }
        int split = 1 + random.nextInt(taxa.size() - 1);
        return "(" + randomNewick(taxa.subList(0, split), random) + "," + randomNewick(taxa.subList(split, taxa.size()), random)
                + "):" + (0.02 + 0.2 * random.nextDouble());
    }

    private static TreeModel randomPhylogeny(int taxa, java.util.Random random) throws Exception {
        java.util.List<String> names = new java.util.ArrayList<>();
        for (int i = 0; i < taxa; i++) {
            names.add("t" + i);
        }
        return new DefaultTreeModel("phylogeny", new NewickImporter(randomNewick(names, random) + ";").importTree(null));
    }

    private static final class SparseProblem {
        final GeneralDataType dataType;
        final String[] codes;
        final TreeModel phylogeny;
        final TreeModel characterTree;
        final Parameter rates;
        final FrequencyModel frequencies;

        SparseProblem(int states, int taxa, java.util.Random random, String[][] ambiguities) throws Exception {
            codes = new String[states];
            for (int i = 0; i < states; i++) {
                codes[i] = "s" + i;
            }
            dataType = new GeneralDataType(codes);
            for (String[] ambiguity : ambiguities) {
                dataType.addAmbiguity(ambiguity[0], Arrays.copyOfRange(ambiguity, 1, ambiguity.length));
            }
            phylogeny = randomPhylogeny(taxa, random);
            characterTree = new DefaultTreeModel("characterTree",
                    new NewickImporter(balancedNewick(codes, 0, states) + ";").importTree(null));
            double[] lambda = new double[states - 1];
            for (int i = 0; i < lambda.length; i++) {
                lambda[i] = 0.3 + 1.7 * random.nextDouble();
            }
            rates = new Parameter.Default(lambda);
            double[] pi = new double[states];
            double total = 0.0;
            for (int i = 0; i < states; i++) {
                pi[i] = 0.5 + random.nextDouble();
                total += pi[i];
            }
            for (int i = 0; i < states; i++) {
                pi[i] /= total;
            }
            frequencies = new FrequencyModel(dataType, new Parameter.Default(pi));
        }

        HierarchicalSubstitutionModel model(String name, HierarchicalSubstitutionModel.Aggregation aggregation) {
            HierarchicalSubstitutionModel model = generalModel(name, characterTree, rates.getParameterValues(), frequencies);
            model.setAggregation(aggregation);
            return model;
        }

        HierarchicalSubstitutionModel shared(String name, HierarchicalSubstitutionModel.Aggregation aggregation) {
            HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel(name,
                    new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, characterTree, frequencies, false));
            model.setAggregation(aggregation);
            return model;
        }

        // 1 to 4 observed codes per column, a fraction of '?' (or of the extra codes) at the tips
        Patterns patterns(int columns, double unknown, int[] extraCodes, double extra, java.util.Random random) {
            Patterns patterns = new Patterns(dataType, phylogeny);
            final int taxa = phylogeny.getTaxonCount();
            for (int k = 0; k < columns; k++) {
                int[] pool = new int[1 + random.nextInt(4)];
                for (int j = 0; j < pool.length; j++) {
                    pool[j] = random.nextInt(codes.length);
                }
                int[] pattern = new int[taxa];
                for (int t = 0; t < taxa; t++) {
                    double u = random.nextDouble();
                    if (u < unknown) {
                        pattern[t] = dataType.getUnknownState();
                    } else if (u < unknown + extra && extraCodes.length > 0) {
                        pattern[t] = extraCodes[random.nextInt(extraCodes.length)];
                    } else {
                        pattern[t] = pool[random.nextInt(pool.length)];
                    }
                }
                patterns.addPattern(pattern);
            }
            return patterns;
        }
    }

    private static BeagleDataLikelihoodDelegate delegate(TreeModel tree, PatternList patterns, SubstitutionModel model,
                                                         GammaSiteRateModel sites, boolean useAmbiguities,
                                                         boolean useHierarchicalRepresentation) {
        return new BeagleDataLikelihoodDelegate(tree, patterns, new HomogeneousBranchModel(model), sites, useAmbiguities,
                false, PartialsRescalingScheme.DEFAULT, false,
                new PreOrderSettings(false, false, false, true, false, false, false, useHierarchicalRepresentation));
    }

    private static final class Paths {
        final BeagleDataLikelihoodDelegate[] delegates = new BeagleDataLikelihoodDelegate[4];
        final TreeDataLikelihood[] likelihoods = new TreeDataLikelihood[4];

        // ON, OFF and AUTO through the hierarchical representation, then the eigen path
        Paths(SparseProblem problem, PatternList patterns, GammaSiteRateModel sites, boolean useAmbiguities) {
            for (int i = 0; i < 4; i++) {
                SubstitutionModel model = problem.shared("hierarchical-" + i, i < 3 ? AGGREGATIONS[i]
                        : HierarchicalSubstitutionModel.Aggregation.AUTO);
                delegates[i] = delegate(problem.phylogeny, patterns, model, sites, useAmbiguities, i < 3);
                likelihoods[i] = new TreeDataLikelihood(delegates[i], problem.phylogeny, new DefaultBranchRateModel());
            }
        }

        double[] evaluate() {
            double[] out = new double[4];
            for (int i = 0; i < 4; i++) {
                out[i] = likelihoods[i].getLogLikelihood();
            }
            return out;
        }

        void store() {
            for (TreeDataLikelihood l : likelihoods) l.storeModelState();
        }

        void restore() {
            for (TreeDataLikelihood l : likelihoods) l.restoreModelState();
        }

        void accept() {
            for (TreeDataLikelihood l : likelihoods) l.acceptModelState();
        }
    }

    private static void assertPathsAgree(String what, double[] values, double expected) {
        final String[] names = {"aggregate=on", "aggregate=off", "aggregate=auto", "eigen path"};
        for (int i = 0; i < values.length; i++) {
            assertEquals(what + " " + names[i], expected, values[i], 1e-10 * Math.abs(expected));
        }
    }

    private static double freshEigen(SparseProblem problem, PatternList patterns, GammaSiteRateModel sites,
                                     boolean useAmbiguities) {
        SubstitutionModel fresh = problem.model("fresh", HierarchicalSubstitutionModel.Aggregation.AUTO);
        return new TreeDataLikelihood(delegate(problem.phylogeny, patterns, fresh, sites, useAmbiguities, false),
                problem.phylogeny, new DefaultBranchRateModel()).getLogLikelihood();
    }

    private void checkSparseChain(String label, SparseProblem problem, PatternList patterns, boolean useAmbiguities,
                                  java.util.Random random) {
        for (int categories : new int[]{1, 4}) {
            final String what = label + " C=" + categories;
            GammaSiteRateModel sites = categories == 1 ? new GammaSiteRateModel("gamma") : gamma(0.5, categories);
            Paths paths = new Paths(problem, patterns, sites, useAmbiguities);
            final double start = freshEigen(problem, patterns, sites, useAmbiguities);
            assertPathsAgree(what + " start", paths.evaluate(), start);

            double[] derivative = null;
            for (int i = 0; i < 4; i++) {
                double[] out = new double[problem.codes.length];
                paths.delegates[i].getRootFrequencyDerivatives(problem.phylogeny.getRoot().getNumber(), out);
                if (derivative == null) {
                    derivative = out;
                    continue;
                }
                for (int x = 0; x < out.length; x++) {
                    assertEquals(what + " root frequency derivative " + x + " of path " + i, derivative[x], out[x],
                            1e-10 * Math.max(1.0, Math.abs(derivative[x])));
                }
            }

            paths.store();
            final double[] lambdaBefore = problem.rates.getParameterValues();
            final int node = random.nextInt(problem.rates.getDimension() - 1);
            final int root = problem.rates.getDimension() - 1;
            problem.rates.setParameterValue(node, lambdaBefore[node] * 1.7);
            problem.rates.setParameterValue(root, lambdaBefore[root] * 0.6);
            final double lambdaMoved = freshEigen(problem, patterns, sites, useAmbiguities);
            assertTrue(what + " the lambda change matters", Math.abs(lambdaMoved - start) > 1e-6);
            assertPathsAgree(what + " lambda change", paths.evaluate(), lambdaMoved);
            paths.restore();
            assertTrue(what + " lambda parameter restored",
                    Arrays.equals(lambdaBefore, problem.rates.getParameterValues()));
            assertPathsAgree(what + " lambda restored", paths.evaluate(), start);

            paths.store();
            swapTips(problem.characterTree, random);
            final double treeMoved = freshEigen(problem, patterns, sites, useAmbiguities);
            assertTrue(what + " the character-tree change matters", Math.abs(treeMoved - start) > 1e-6);
            assertPathsAgree(what + " character-tree change", paths.evaluate(), treeMoved);
            paths.accept();
            for (TreeDataLikelihood l : paths.likelihoods) l.makeDirty();
            assertPathsAgree(what + " full recompute", paths.evaluate(), treeMoved);
        }
    }

    public void testSparseGeneralStatesAggregation() throws Exception {
        if (skip("testSparseGeneralStatesAggregation")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261104);
        SparseProblem problem = new SparseProblem(256, 40, random, new String[0][]);
        Patterns patterns = problem.patterns(200, 0.05, new int[0], 0.0, random);
        checkSparseChain("S=256 states", problem, patterns, false, random);
    }

    public void testSparseAmbiguityCodesAggregation() throws Exception {
        if (skip("testSparseAmbiguityCodesAggregation")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261105);
        SparseProblem problem = new SparseProblem(256, 40, random,
                new String[][]{{"near", "s0", "s1"}, {"far", "s3", "s200", "s201"}});
        int[] extra = {problem.dataType.getState("near"), problem.dataType.getState("far")};
        Patterns patterns = problem.patterns(200, 0.05, extra, 0.05, random);
        checkSparseChain("S=256 ambiguity codes, partials", problem, patterns, true, random);
        checkSparseChain("S=256 ambiguity codes, states", problem, patterns, false, random);
    }

    public void testSparseUncertainColumnsAggregation() throws Exception {
        if (skip("testSparseUncertainColumnsAggregation")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261106);
        SparseProblem problem = new SparseProblem(64, 24, random, new String[0][]);
        final int S = problem.codes.length, taxa = problem.phylogeny.getTaxonCount();
        UncertainSiteList patterns = new UncertainSiteList(problem.dataType, problem.phylogeny);
        for (int k = 0; k < 60; k++) {
            double[][] column = new double[taxa][S];
            if (k == 0) {
                for (double[] tip : column) Arrays.fill(tip, 0.5);
            } else {
                int a = random.nextInt(S), b = random.nextInt(S);
                for (double[] tip : column) {
                    double u = random.nextDouble();
                    if (u < 0.05) {
                        Arrays.fill(tip, 1.0);
                    } else if (u < 0.15) {
                        tip[a] = 0.3;
                        tip[b] = 0.7;
                    } else {
                        tip[u < 0.6 ? a : b] = 1.0;
                    }
                }
            }
            patterns.addPattern(column);
        }
        checkSparseChain("S=64 uncertain columns", problem, patterns, true, random);
    }

    // Blocked aggregation at 4 states: BEAST orders the patterns by aggregated block (HierarchicalRateAggregatedPattern)
    // and the plugin runs each block with a reduced kernel. The log-likelihood must not change.
    public void testBlockSortedPatternsThroughBeagle() throws Exception {
        if (skip("testBlockSortedPatternsThroughBeagle")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261119);
        TreeModel phylogeny = randomPhylogeny(40, random);
        final int[][] pools = {{0, 2}, {1, 3}, {0, 1, 2, 3}};
        Patterns patterns = new Patterns(Nucleotides.INSTANCE, phylogeny);
        for (int k = 0; k < 900; k++) {
            int[] pool = pools[random.nextInt(3)];
            int[] pattern = new int[phylogeny.getTaxonCount()];
            for (int t = 0; t < pattern.length; t++) {
                pattern[t] = random.nextDouble() < 0.05 ? Nucleotides.INSTANCE.getUnknownState()
                        : pool[random.nextInt(pool.length)];
            }
            patterns.addPattern(pattern, 1.0 + random.nextInt(3));
        }
        final double[] pi = {0.3, 0.2, 0.25, 0.25};
        final String tn93 = "((A:1,G:1):1,(C:1,T:1):1);";
        HierarchicalSubstitutionModel on = hierarchical(tn93, pi, 1.2, 2.3, 0.9);
        on.setAggregation(HierarchicalSubstitutionModel.Aggregation.ON);
        HierarchicalSubstitutionModel off = hierarchical(tn93, pi, 1.2, 2.3, 0.9);
        off.setAggregation(HierarchicalSubstitutionModel.Aggregation.OFF);

        HierarchicalRateAggregatedPattern sorted = new HierarchicalRateAggregatedPattern(patterns, on.getBeagleModelTriples());
        int[] sizes = sorted.getBlockSizes();
        assertEquals(patterns.getPatternCount(), sorted.getPatternCount());
        assertEquals(patterns.getPatternCount(), sizes[0] + sizes[1] + sizes[2]);
        assertTrue("both reduced blocks are present", sizes[0] > 0 && sizes[2] > 0);
        java.util.Map<String, Double> before = new java.util.HashMap<>(), after = new java.util.HashMap<>();
        for (int k = 0; k < patterns.getPatternCount(); k++) {
            before.merge(Arrays.toString(patterns.getPattern(k)), patterns.getPatternWeight(k), Double::sum);
            after.merge(Arrays.toString(sorted.getPattern(k)), sorted.getPatternWeight(k), Double::sum);
        }
        assertEquals("same patterns and weights", before, after);
        int group = 0;
        for (int k = 0; k < sorted.getPatternCount(); k++) {
            int seen = 0;
            for (int state : sorted.getPattern(k)) {
                if (state < S) {
                    seen |= 1 << state;
                }
            }
            int g = seen != 0 && (seen & ~5) == 0 ? 0 : seen != 0 && (seen & ~10) == 0 ? 2 : 1;
            assertTrue("blocks are contiguous and ordered at pattern " + k, g >= group);
            group = g;
        }

        assertSame(patterns, HierarchicalSubstitutionModelParser.aggregatedPatterns(patterns, new HomogeneousBranchModel(off)));
        assertTrue(HierarchicalSubstitutionModelParser.aggregatedPatterns(patterns, new HomogeneousBranchModel(on))
                instanceof HierarchicalRateAggregatedPattern);

        GammaSiteRateModel sites = gamma(0.5, 4);
        TreeDataLikelihood reference = new TreeDataLikelihood(delegate(phylogeny, patterns, off, sites, false, true),
                phylogeny, new DefaultBranchRateModel());
        TreeDataLikelihood blocked = new TreeDataLikelihood(delegate(phylogeny, sorted, on, sites, false, true),
                phylogeny, new DefaultBranchRateModel());
        TreeDataLikelihood unsorted = new TreeDataLikelihood(delegate(phylogeny, patterns, on, sites, false, true),
                phylogeny, new DefaultBranchRateModel());
        for (int round = 0; round < 3; round++) {
            if (round > 0) {
                NodeRef node = phylogeny.getNode(phylogeny.getExternalNodeCount() + random.nextInt(
                        phylogeny.getInternalNodeCount() - 1));
                double height = phylogeny.getNodeHeight(node);
                double low = 0.0;
                for (int c = 0; c < phylogeny.getChildCount(node); c++) {
                    low = Math.max(low, phylogeny.getNodeHeight(phylogeny.getChild(node, c)));
                }
                phylogeny.setNodeHeight(node, low + 0.9 * (height - low));
            }
            double expected = reference.getLogLikelihood();
            assertEquals("sorted, aggregation on, round " + round, expected, blocked.getLogLikelihood(),
                    1e-10 * Math.abs(expected));
            assertEquals("unsorted, aggregation on, round " + round, expected, unsorted.getLogLikelihood(),
                    1e-10 * Math.abs(expected));
        }
    }

    // Blocked aggregation at 4 states under a sampled character tree: NNI, keep and rate moves on block-sized data,
    // with the patterns sorted for the starting tree and unsorted. Aggregation on must match off at every step.
    public void testBlockedCharacterTreeChainThroughBeagle() throws Exception {
        if (skip("testBlockedCharacterTreeChainThroughBeagle")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261120);
        TreeModel phylogeny = randomPhylogeny(40, random);
        final int[][] pools = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}, {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3},
                {0, 1, 2, 3}};
        Patterns patterns = new Patterns(Nucleotides.INSTANCE, phylogeny);
        for (int k = 0; k < 2400; k++) {
            int[] pool = pools[random.nextInt(pools.length)];
            int[] pattern = new int[phylogeny.getTaxonCount()];
            for (int t = 0; t < pattern.length; t++) {
                pattern[t] = random.nextDouble() < 0.05 ? Nucleotides.INSTANCE.getUnknownState()
                        : pool[random.nextInt(pool.length)];
            }
            patterns.addPattern(pattern, 1.0 + random.nextInt(3));
        }
        final FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(new double[]{0.3, 0.2, 0.25, 0.25}));
        Parameter keep = new Parameter.Default(1);
        dr.evomodel.substmodel.CharacterTreeModel tree = new dr.evomodel.substmodel.CharacterTreeModel("characterTree",
                "((A,G),(C,T));", Nucleotides.INSTANCE, true, keep);
        Parameter rates = new Parameter.Default(new double[]{1.4, 0.7, 1.1});
        final HierarchicalSubstitutionModel.Aggregation[] aggregations = {HierarchicalSubstitutionModel.Aggregation.ON,
                HierarchicalSubstitutionModel.Aggregation.ON, HierarchicalSubstitutionModel.Aggregation.OFF};
        HierarchicalSubstitutionModel[] models = new HierarchicalSubstitutionModel[3];
        for (int i = 0; i < 3; i++) {
            models[i] = new HierarchicalSubstitutionModel("hierarchical-" + i,
                    new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, tree.getContraction(), frequencies, true));
            models[i].setAggregation(aggregations[i]);
        }
        HierarchicalRateAggregatedPattern sorted = new HierarchicalRateAggregatedPattern(patterns, models[0].getBeagleModelTriples());
        int[] sizes = sorted.getBlockSizes();
        assertTrue("both reduced blocks hold at least 150 patterns: " + Arrays.toString(sizes),
                sizes[0] >= 150 && sizes[2] >= 150);
        final PatternList[] lists = {sorted, patterns, patterns};
        GammaSiteRateModel sites = gamma(0.6, 4);
        TreeDataLikelihood[] paths = new TreeDataLikelihood[3];
        for (int i = 0; i < 3; i++) {
            paths[i] = new TreeDataLikelihood(delegate(phylogeny, lists[i], models[i], sites, false, true), phylogeny,
                    new DefaultBranchRateModel());
        }
        final String[] names = {"aggregate=on, sorted", "aggregate=on, unsorted", "aggregate=off"};
        final dr.evomodel.tree.DefaultTreeModel skeleton = tree.getSkeleton();
        final int root = skeleton.getRoot().getNumber();
        java.util.Set<String> visited = new java.util.HashSet<>();
        int[] kinds = new int[4];
        final String start = tree.getContractedNewick();
        int left = 0, returned = 0, atStart = 0;
        double current = paths[2].getLogLikelihood();
        for (int i = 0; i < 2; i++) {
            assertEquals(names[i] + " start", current, paths[i].getLogLikelihood(), 1e-10 * Math.abs(current));
        }
        for (int step = 0; step < 300; step++) {
            final String before = tree.getContractedNewick();
            for (TreeDataLikelihood l : paths) l.storeModelState();
            tree.storeModelState();
            final int kind = random.nextInt(4);
            kinds[kind]++;
            int slot;
            do {
                slot = random.nextInt(S - 1);
            } while (S + slot == root);
            if (kind == 0) {
                NodeRef u = skeleton.getNode(S + slot);
                NodeRef p = skeleton.getParent(u);
                NodeRef w = skeleton.getChild(p, 0) == u ? skeleton.getChild(p, 1) : skeleton.getChild(p, 0);
                NodeRef c = skeleton.getChild(u, random.nextInt(2));
                skeleton.beginTreeEdit();
                skeleton.removeChild(u, c);
                skeleton.removeChild(p, w);
                skeleton.addChild(u, w);
                skeleton.addChild(p, c);
                tree.resetUnitHeightsQuietly();
                skeleton.endTreeEdit();
            } else if (kind == 1) {
                if (keep.getParameterValue(slot) == 0.0) {
                    rates.setParameterValue(slot, 0.05 + random.nextDouble());
                    keep.setParameterValue(slot, 1.0);
                } else {
                    keep.setParameterValue(slot, 0.0);
                    rates.setParameterValue(slot, 0.05 + random.nextDouble());
                }
            } else if (kind == 2) {
                final int r = random.nextInt(S - 1);
                rates.setParameterValue(r, rates.getParameterValue(r) * Math.exp(0.5 * (random.nextDouble() - 0.5)));
            } else {
                NodeRef node = phylogeny.getNode(phylogeny.getExternalNodeCount() + random.nextInt(
                        phylogeny.getInternalNodeCount() - 1));
                double height = phylogeny.getNodeHeight(node);
                double low = 0.0;
                for (int c = 0; c < phylogeny.getChildCount(node); c++) {
                    low = Math.max(low, phylogeny.getNodeHeight(phylogeny.getChild(node, c)));
                }
                phylogeny.setNodeHeight(node, low + (0.8 + 0.15 * random.nextDouble()) * (height - low));
            }
            final String now = tree.getContractedNewick();
            visited.add(now);
            final double expected = paths[2].getLogLikelihood();
            HierarchicalSubstitutionModel fresh = sampledModel(tree, new Parameter.Default(rates.getParameterValues()),
                    new Parameter.Default(keep.getParameterValues()), frequencies);
            fresh.setAggregation(HierarchicalSubstitutionModel.Aggregation.OFF);
            final double reference = new TreeDataLikelihood(delegate(phylogeny, patterns, fresh, sites, false, true),
                    phylogeny, new DefaultBranchRateModel()).getLogLikelihood();
            assertEquals("step " + step + " aggregate=off against a fresh model", reference, expected,
                    1e-10 * Math.abs(reference));
            for (int i = 0; i < 2; i++) {
                assertEquals("step " + step + " " + names[i], expected, paths[i].getLogLikelihood(),
                        1e-10 * Math.abs(expected));
            }
            if (now.equals(start) || (before.equals(start) ? random.nextDouble() < 0.4 : kind >= 2 && random.nextBoolean())) {
                for (TreeDataLikelihood l : paths) l.acceptModelState();
                tree.acceptModelState();
                current = expected;
                if (now.equals(start)) {
                    atStart++;
                    if (!before.equals(start)) returned++;
                } else if (before.equals(start)) {
                    left++;
                }
            } else {
                for (TreeDataLikelihood l : paths) l.restoreModelState();
                tree.restoreModelState();
                for (int i = 0; i < 3; i++) {
                    assertEquals("step " + step + " " + names[i] + " restored", current, paths[i].getLogLikelihood(),
                            1e-10 * Math.abs(current));
                }
            }
        }
        assertTrue("character trees visited " + visited, visited.size() >= 8);
        assertTrue("left the starting tree " + left + " times, returned " + returned + ", accepted there " + atStart,
                left >= 10 && returned >= 10 && atStart >= 50);
        for (int k : kinds) assertTrue("every move kind ran " + Arrays.toString(kinds), k > 20);
    }

    public void testAggregationSummary() throws Exception {
        String[] codes = new String[8];
        for (int i = 0; i < 8; i++) {
            codes[i] = "s" + i;
        }
        GeneralDataType dataType = new GeneralDataType(codes);
        TreeModel characterTree = new DefaultTreeModel("characterTree",
                new NewickImporter(balancedNewick(codes, 0, 8) + ";").importTree(null));
        FrequencyModel frequencies = new FrequencyModel(dataType, new Parameter.Default(8, 0.125));
        HierarchicalSubstitutionModel model = generalModel("summary", characterTree, new double[]{1, 1, 1, 1, 1, 1, 1},
                frequencies);
        boolean[][] observed = new boolean[5][8];
        observed[1][0] = true;
        observed[2][0] = observed[2][1] = true;
        observed[3][0] = observed[3][7] = true;
        observed[4][0] = observed[4][1] = true;
        HierarchicalSubstitutionModel.AggregationSummary summary = HierarchicalSubstitutionModel.summarizeAggregation(
                model.getBeagleModelTriples(), 8, 1, observed);
        assertEquals(4, summary.classCount);
        assertEquals(1 + 4 + 4 + 6 + 4, summary.reducedStateSum);
        assertEquals(6, summary.maxReducedStates);
        assertEquals(19.0 / 40.0, summary.getRatio(), 1e-15);
        assertFalse(summary.automaticRuleAggregates());
        assertTrue(summary.aggregates(HierarchicalSubstitutionModel.Aggregation.ON));

        boolean[] marks = new boolean[4];
        HierarchicalSubstitutionModel.markObservedStates(new double[]{0.5, 0.5, 0.5, 0.5}, 0, 4, marks);
        assertFalse(marks[0] || marks[1] || marks[2] || marks[3]);
        HierarchicalSubstitutionModel.markObservedStates(new double[]{0.0, 0.3, 0.0, 0.7}, 0, 4, marks);
        assertTrue(!marks[0] && marks[1] && !marks[2] && marks[3]);
        Arrays.fill(marks, false);
        HierarchicalSubstitutionModel.markObservedStates(new double[]{0.1, 0.3, 0.2, 0.7}, 0, 4, marks);
        assertTrue(marks[0] && marks[1] && marks[2] && marks[3]);
    }

    public void testSeparateRootFrequenciesRefused() throws Exception {
        if (skip("testSeparateRootFrequenciesRefused")) {
            return;
        }
        HierarchicalSubstitutionModel model = hierarchical("((A:1,G:1):1,(C:1,T:1):1);",
                new double[]{0.1, 0.2, 0.3, 0.4}, 1.0, 2.0, 3.0);
        FrequencyModel root = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25}));
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        try {
            new BeagleDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(model, root),
                    new GammaSiteRateModel("gamma"), false, false, PartialsRescalingScheme.DEFAULT, false,
                    new PreOrderSettings(false, false, false, true, false, false, false, true));
            fail("a separate root frequency model must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("root frequencies"));
        }
        TreeDataLikelihood same = new TreeDataLikelihood(new BeagleDataLikelihoodDelegate(treeModel, patterns,
                new HomogeneousBranchModel(model, model.getFrequencyModel()), new GammaSiteRateModel("gamma"), false,
                false, PartialsRescalingScheme.DEFAULT, false,
                new PreOrderSettings(false, false, false, true, false, false, false, true)), treeModel,
                new DefaultBranchRateModel());
        assertEquals(likelihood(model, new GammaSiteRateModel("gamma"), false).getLogLikelihood(),
                same.getLogLikelihood(), 1e-10 * Math.abs(same.getLogLikelihood()));
    }

    public void testAggregateFlagsReachBeagle() throws Exception {
        if (skip("testAggregateFlagsReachBeagle")) {
            return;
        }
        final String override = System.getenv("BEAGLE_HRATE_AGGREGATE");
        java.util.Random random = new java.util.Random(20261107);
        SparseProblem problem = new SparseProblem(32, 8, random, new String[0][]);
        Patterns patterns = problem.patterns(10, 0.05, new int[0], 0.0, random);
        java.util.logging.Logger logger = java.util.logging.Logger.getLogger("dr.evomodel");
        final StringBuilder messages = new StringBuilder();
        java.util.logging.Handler handler = new java.util.logging.Handler() {
            public void publish(java.util.logging.LogRecord record) {
                messages.append(record.getMessage()).append('\n');
            }
            public void flush() {
            }
            public void close() {
            }
        };
        logger.addHandler(handler);
        try {
            for (HierarchicalSubstitutionModel.Aggregation aggregation : AGGREGATIONS) {
                messages.setLength(0);
                delegate(problem.phylogeny, patterns, problem.model("flags", aggregation), new GammaSiteRateModel("gamma"),
                        false, true);
                final String log = messages.toString();
                final String flags = log.substring(log.indexOf("with instance flags:"));
                final String line = flags.substring(0, flags.indexOf('\n'));
                HierarchicalSubstitutionModel.Aggregation effective = "1".equals(override)
                        ? HierarchicalSubstitutionModel.Aggregation.ON : "0".equals(override)
                        ? HierarchicalSubstitutionModel.Aggregation.OFF : aggregation;
                assertEquals(aggregation + ": " + line, effective == HierarchicalSubstitutionModel.Aggregation.ON,
                        line.contains(BeagleFlag.HIERARCHICAL_AGGREGATE_ON.name()));
                assertEquals(aggregation + ": " + line, effective == HierarchicalSubstitutionModel.Aggregation.OFF,
                        line.contains(BeagleFlag.HIERARCHICAL_AGGREGATE_OFF.name()));
                assertTrue(log, log.contains("aggregate=" + aggregation.name().toLowerCase() + ":"));
            }
        } finally {
            logger.removeHandler(handler);
        }
    }

    public void testCodonCharacterTree() throws Exception {
        dr.evolution.datatype.Codons codons = dr.evolution.datatype.Codons.UNIVERSAL;
        TreeModel characterTree = AbstractHierarchicalRateProvider.defaultCharacterTree(codons);
        assertEquals(61, characterTree.getExternalNodeCount());
        for (int state = 0; state < codons.getStateCount(); state++) {
            final String code = AbstractHierarchicalRateProvider.codeOfState(codons, state);
            assertEquals(3, code.length());
            assertEquals(state, AbstractHierarchicalRateProvider.stateOfCode(codons, code));
        }
        assertEquals(-1, AbstractHierarchicalRateProvider.stateOfCode(codons, "TAA"));
        assertEquals(2, AbstractHierarchicalRateProvider.stateOfCode(Nucleotides.INSTANCE, "G"));
        double[] pi = new double[61];
        Arrays.fill(pi, 1.0 / 61);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("codon",
                new EqualRatesProvider(characterTree, new FrequencyModel(codons, new Parameter.Default(pi)), true));
        double[] matrix = new double[61 * 61];
        model.getTransitionProbabilities(0.1, matrix);
        for (int x = 0; x < 61; x++) {
            double row = 0.0;
            for (int y = 0; y < 61; y++) {
                row += matrix[x * 61 + y];
            }
            assertEquals(1.0, row, 1e-12);
        }
    }

    public void testMultifurcatingTreesThroughBeagle() throws Exception {
        if (skip("testMultifurcatingTreesThroughBeagle")) {
            return;
        }
        FrequencyModel frequencies = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        dr.evomodel.substmodel.CharacterTreeContraction.Built star = dr.evomodel.substmodel.CharacterTreeContraction.build(
                "characterTree", "(A,C,G,T);", Nucleotides.INSTANCE, true, null);
        HierarchicalSubstitutionModel f81 = new HierarchicalSubstitutionModel("hierarchical", new EqualRatesProvider(
                new dr.evomodel.substmodel.CharacterTreeContraction(star.skeleton, Nucleotides.INSTANCE, star.collapsed, null),
                frequencies, true));
        HKY hky1 = new HKY(new Parameter.Default(1, 1.0), frequencies);
        final double expected = likelihood(hky1, gamma(0.5, 4), false).getLogLikelihood();
        for (HierarchicalSubstitutionModel.Aggregation aggregation : AGGREGATIONS) {
            f81.setAggregation(aggregation);
            assertEquals("star = F81, aggregate " + aggregation, expected,
                    likelihood(f81, gamma(0.5, 4), true).getLogLikelihood(), 1e-10 * Math.abs(expected));
        }
        assertEquals("star on the eigen path", expected, likelihood(f81, gamma(0.5, 4), false).getLogLikelihood(),
                1e-10 * Math.abs(expected));

        dr.evomodel.substmodel.CharacterTreeContraction.Built oneClade = dr.evomodel.substmodel.CharacterTreeContraction.build(
                "characterTree", "((A,G),C,T);", Nucleotides.INSTANCE, true, null);
        Parameter rates = new Parameter.Default(oneClade.skeleton.getInternalNodeCount(), 0.0);
        final double rootRate = 1.1, agRate = 0.8;
        for (int i = 0; i < rates.getDimension(); i++) {
            String c = clade(oneClade.skeleton, oneClade.skeleton.getNode(S + i));
            rates.setParameterValue(i, c.equals("AG") ? agRate : c.equals("ACGT") ? rootRate : 3.3);
        }
        HierarchicalSubstitutionModel hierarchical = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates,
                        new dr.evomodel.substmodel.CharacterTreeContraction(oneClade.skeleton, Nucleotides.INSTANCE,
                                oneClade.collapsed, null), frequencies, true));
        double[] pi = frequencies.getFrequencies();
        TN93 tn93 = new TN93(new Parameter.Default(1, (rootRate + agRate / (pi[0] + pi[2])) / rootRate),
                new Parameter.Default(1, 1.0), frequencies);
        final double tn93Expected = likelihood(tn93, gamma(0.5, 4), false).getLogLikelihood();
        assertEquals("((A,G),C,T) = TN93(kappa1, 1)", tn93Expected,
                likelihood(hierarchical, gamma(0.5, 4), true).getLogLikelihood(), 1e-10 * Math.abs(tn93Expected));
    }

    private static String randomMultifurcatingNewick(java.util.List<String> codes, java.util.Random random) {
        java.util.List<String> pieces = new java.util.ArrayList<>(codes);
        while (pieces.size() > 4) {
            final int k = 2 + random.nextInt(3);
            java.util.List<String> chosen = new java.util.ArrayList<>();
            for (int i = 0; i < k; i++) chosen.add(pieces.remove(random.nextInt(pieces.size())));
            pieces.add(random.nextInt(pieces.size() + 1), "(" + String.join(",", chosen) + ")");
        }
        return "(" + String.join(",", pieces) + ");";
    }

    private static void exchangeWithUncle(TreeModel tree, NodeRef u, int which) {
        NodeRef p = tree.getParent(u);
        NodeRef w = tree.getChild(p, 0) == u ? tree.getChild(p, 1) : tree.getChild(p, 0);
        NodeRef c = tree.getChild(u, which);
        tree.beginTreeEdit();
        tree.removeChild(u, c);
        tree.removeChild(p, w);
        tree.addChild(u, w);
        tree.addChild(p, c);
        tree.endTreeEdit();
    }

    public void testNeutralCharacterTreeMovesThroughBeagle() throws Exception {
        if (skip("testNeutralCharacterTreeMovesThroughBeagle")) {
            return;
        }
        java.util.Random random = new java.util.Random(20261114);
        for (int states : new int[]{8, 12, 20}) {
            SparseProblem problem = new SparseProblem(states, 12, random, new String[0][]);
            Patterns patterns = problem.patterns(120, 0.1, new int[0], 0.0, random);
            String newick = randomMultifurcatingNewick(Arrays.asList(problem.codes), random);
            dr.evomodel.substmodel.CharacterTreeContraction.Built built = dr.evomodel.substmodel.CharacterTreeContraction.build(
                    "characterTree", newick, problem.dataType, false, null);
            final TreeModel skeleton = built.skeleton;
            assertTrue(newick, built.isMultifurcating());
            double[] raw = new double[skeleton.getInternalNodeCount()];
            for (int i = 0; i < raw.length; i++) raw[i] = 0.3 + random.nextDouble();
            Parameter rates = new Parameter.Default(raw.clone());
            dr.evomodel.substmodel.CharacterTreeContraction contraction = new dr.evomodel.substmodel.CharacterTreeContraction(
                    skeleton, problem.dataType, built.collapsed, null);
            GammaSiteRateModel sites = gamma(0.7, 3);
            TreeDataLikelihood[] paths = new TreeDataLikelihood[4];
            for (int i = 0; i < 4; i++) {
                HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical-" + i,
                        new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, contraction, problem.frequencies, true));
                model.setAggregation(i < 3 ? AGGREGATIONS[i] : HierarchicalSubstitutionModel.Aggregation.AUTO);
                paths[i] = new TreeDataLikelihood(delegate(problem.phylogeny, patterns, model, sites, false, i < 3),
                        problem.phylogeny, new DefaultBranchRateModel());
            }
            double[] first = new double[4];
            for (int i = 0; i < 4; i++) first[i] = paths[i].getLogLikelihood();
            assertPathsAgree("S=" + states + " start", first, first[3]);
            final String contracted = contraction.contractedNewick(null);

            int neutral = 0;
            for (int m = 0; m < 40 && neutral < 6; m++) {
                NodeRef u = skeleton.getNode(states + random.nextInt(states - 1));
                if (skeleton.isRoot(u)) continue;
                boolean isNeutral = contraction.isCollapsed(u.getNumber());
                for (TreeDataLikelihood l : paths) l.storeModelState();
                exchangeWithUncle(skeleton, u, random.nextInt(2));
                double[] now = new double[4];
                for (int i = 0; i < 4; i++) now[i] = paths[i].getLogLikelihood();
                HierarchicalSubstitutionModel fresh = new HierarchicalSubstitutionModel("fresh",
                        new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(raw.clone()),
                                new dr.evomodel.substmodel.CharacterTreeContraction(skeleton, problem.dataType, built.collapsed, null),
                                problem.frequencies, true));
                fresh.setAggregation(HierarchicalSubstitutionModel.Aggregation.OFF);
                final double freshHierarchical = new TreeDataLikelihood(delegate(problem.phylogeny, patterns, fresh, sites, false, true),
                        problem.phylogeny, new DefaultBranchRateModel()).getLogLikelihood();
                if (isNeutral) {
                    neutral++;
                    assertEquals(contracted, contraction.contractedNewick(null));
                    for (int i = 0; i < 4; i++) {
                        assertEquals("S=" + states + " neutral move, path " + i, first[i], now[i], 0.0);
                    }
                    assertEquals("S=" + states + " fresh model after a neutral move", first[1], freshHierarchical, 0.0);
                    for (TreeDataLikelihood l : paths) l.acceptModelState();
                } else {
                    assertPathsAgree("S=" + states + " after a move", now, freshHierarchical);
                    for (TreeDataLikelihood l : paths) l.restoreModelState();
                    for (int i = 0; i < 4; i++) {
                        assertEquals("S=" + states + " restored, path " + i, first[i], paths[i].getLogLikelihood(),
                                i < 3 ? 0.0 : 1e-10 * Math.abs(first[i]));
                    }
                }
            }
            assertTrue("S=" + states + " neutral moves tried", neutral > 0);
        }
    }

    private static HierarchicalSubstitutionModel sampledModel(dr.evomodel.substmodel.CharacterTreeModel tree, Parameter rates,
                                                              Parameter keep, FrequencyModel frequencies) {
        return new HierarchicalSubstitutionModel("hierarchical", new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(
                rates, new dr.evomodel.substmodel.CharacterTreeContraction(tree.getSkeleton(), tree.getDataType(), null, keep),
                frequencies, true));
    }

    public void testSampledCharacterTreeChainThroughBeagle() throws Exception {
        if (skip("testSampledCharacterTreeChainThroughBeagle")) {
            return;
        }
        java.util.Random random = new java.util.Random(1114);
        for (int states : new int[]{4, 8}) {
            final TreeModel phylogeny;
            final PatternList patterns;
            final FrequencyModel frequencies;
            final dr.evolution.datatype.DataType type;
            final String newick;
            if (states == 4) {
                phylogeny = treeModel;
                patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
                frequencies = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
                type = Nucleotides.INSTANCE;
                newick = "((A,G),(C,T));";
            } else {
                SparseProblem problem = new SparseProblem(states, 10, random, new String[0][]);
                phylogeny = problem.phylogeny;
                patterns = problem.patterns(100, 0.1, new int[0], 0.0, random);
                frequencies = problem.frequencies;
                type = problem.dataType;
                newick = "((s0,s1,s2),(s3,(s4,s5)),s6,s7);";
            }
            Parameter keep = new Parameter.Default(1);
            dr.evomodel.substmodel.CharacterTreeModel tree = new dr.evomodel.substmodel.CharacterTreeModel("characterTree",
                    newick, type, true, keep);
            Parameter rates = new Parameter.Default(states - 1, 1.0);
            for (int i = 0; i < states - 1; i++) rates.setParameterValue(i, 0.3 + random.nextDouble());
            HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                    new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, tree.getContraction(), frequencies, true));
            GammaSiteRateModel sites = gamma(0.8, 3);
            TreeDataLikelihood hierarchical = new TreeDataLikelihood(delegate(phylogeny, patterns, model, sites, false, true),
                    phylogeny, new DefaultBranchRateModel());
            TreeDataLikelihood eigen = new TreeDataLikelihood(delegate(phylogeny, patterns, model, sites, false, false),
                    phylogeny, new DefaultBranchRateModel());
            final dr.evomodel.tree.DefaultTreeModel skeleton = tree.getSkeleton();
            final int root = skeleton.getRoot().getNumber();
            int[] kinds = new int[3];
            java.util.Set<String> visited = new java.util.HashSet<>();
            double current = hierarchical.getLogLikelihood();
            assertEquals(current, eigen.getLogLikelihood(), 1e-10 * Math.abs(current));
            for (int step = 0; step < 150; step++) {
                hierarchical.storeModelState();
                eigen.storeModelState();
                tree.storeModelState();
                final int kind = random.nextInt(3);
                kinds[kind]++;
                int slot;
                do {
                    slot = random.nextInt(states - 1);
                } while (states + slot == root);
                if (kind == 0) {
                    NodeRef u = skeleton.getNode(states + slot);
                    NodeRef p = skeleton.getParent(u);
                    NodeRef w = skeleton.getChild(p, 0) == u ? skeleton.getChild(p, 1) : skeleton.getChild(p, 0);
                    NodeRef c = skeleton.getChild(u, random.nextInt(2));
                    skeleton.beginTreeEdit();
                    skeleton.removeChild(u, c);
                    skeleton.removeChild(p, w);
                    skeleton.addChild(u, w);
                    skeleton.addChild(p, c);
                    tree.resetUnitHeightsQuietly();
                    skeleton.endTreeEdit();
                } else if (kind == 1) {
                    if (keep.getParameterValue(slot) == 0.0) {
                        rates.setParameterValue(slot, 0.05 + random.nextDouble());
                        keep.setParameterValue(slot, 1.0);
                    } else {
                        keep.setParameterValue(slot, 0.0);
                        rates.setParameterValue(slot, 0.05 + random.nextDouble());
                    }
                } else {
                    rates.setParameterValue(slot, rates.getParameterValue(slot) * Math.exp(0.5 * (random.nextDouble() - 0.5)));
                }
                visited.add(tree.getContractedNewick());
                final double h = hierarchical.getLogLikelihood();
                final double e = eigen.getLogLikelihood();
                Parameter keepCopy = new Parameter.Default(keep.getParameterValues());
                HierarchicalSubstitutionModel fresh = sampledModel(tree, new Parameter.Default(rates.getParameterValues()), keepCopy,
                        frequencies);
                final double expected = new TreeDataLikelihood(delegate(phylogeny, patterns, fresh, sites, false, false),
                        phylogeny, new DefaultBranchRateModel()).getLogLikelihood();
                assertEquals("S=" + states + " step " + step + " hierarchical", expected, h, 1e-10 * Math.abs(expected));
                assertEquals("S=" + states + " step " + step + " eigen", expected, e, 1e-10 * Math.abs(expected));
                if (random.nextBoolean()) {
                    hierarchical.acceptModelState();
                    eigen.acceptModelState();
                    tree.acceptModelState();
                    current = h;
                } else {
                    hierarchical.restoreModelState();
                    eigen.restoreModelState();
                    tree.restoreModelState();
                    assertEquals("S=" + states + " step " + step + " restored", current, hierarchical.getLogLikelihood(),
                            1e-10 * Math.abs(current));
                    assertEquals("S=" + states + " step " + step + " eigen restored", current, eigen.getLogLikelihood(),
                            1e-10 * Math.abs(current));
                }
            }
            assertTrue("S=" + states + " trees visited " + visited.size(), visited.size() > (states == 4 ? 10 : 20));
        }
    }
}
