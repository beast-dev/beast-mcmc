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
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BeagleFunctionality;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.HKYParser;
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

    private static HierarchicalSubstitutionModel hierarchicalHKY(double[] pi, double kappa) throws Exception {
        final int[] partner = {2, 3, 0, 1};
        double norm = 0.0;
        for (int i = 0; i < S; i++) {
            for (int j = 0; j < S; j++) {
                if (i != j) {
                    norm += pi[i] * pi[j] * (j == partner[i] ? kappa : 1.0);
                }
            }
        }
        final double beta = 1.0 / norm;
        return hierarchical("((A:1,G:1):1,(C:1,T:1):1);", pi, beta,
                (kappa - 1.0) * beta * (pi[0] + pi[2]), (kappa - 1.0) * beta * (pi[1] + pi[3]));
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
        HierarchicalSubstitutionModel model = hierarchicalHKY(pi, 8.0);
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0),
                new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(pi)));

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
}
