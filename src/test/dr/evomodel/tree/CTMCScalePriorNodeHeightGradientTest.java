/*
 * CTMCScalePriorNodeHeightGradientTest.java
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

package test.dr.evomodel.tree;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.CTMCScalePrior;
import dr.evomodel.tree.CTMCScalePriorNodeHeightGradient;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Node-height gradient of the CTMC rate reference prior, which depends on the heights through the tree length
 * (plan 3, F3), against central finite differences of the prior.
 */
public class CTMCScalePriorNodeHeightGradientTest extends TestCase {

    // heterochronous
    private static final String NEWICK = "((A:0.12,B:0.32):0.31,(C:0.25,(D:0.07,E:0.17):0.18):0.18);";

    public CTMCScalePriorNodeHeightGradientTest(String name) {
        super(name);
    }

    public void testWholeTree() throws Exception {
        TreeModel tree = tree();
        check(tree, new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, false, null, false), "plain");
    }

    public void testSeveralRates() throws Exception {
        TreeModel tree = tree();
        check(tree, new CTMCScalePrior("prior", new Parameter.Default(new double[]{0.8, 2.5}), tree, null, false,
                null, false), "two rates");
    }

    public void testReciprocal() throws Exception {
        TreeModel tree = tree();
        check(tree, new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, true, null, false),
                "reciprocal");
    }

    public void testSubstitutionModelScaling() throws Exception {
        TreeModel tree = tree();
        check(tree, new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, false, hky(), false),
                "substitution model");
        check(tree, new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, false, hky(), true), "trial");
    }

    public void testSubTree() throws Exception {
        TreeModel tree = tree();
        for (String[] names : new String[][]{{"C", "D", "E"}, {"A", "D", "E"}, {"D", "E"}}) {
            Taxa taxa = new Taxa();
            for (String name : names) {
                taxa.addTaxon(tree.getTaxon(tree.getTaxonIndex(name)));
            }
            check(tree, new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, taxa, false, null, false),
                    "subtree " + Arrays.toString(names));
        }
    }

    public void testTwoTreesWithProportion() throws Exception {
        TreeModel first = tree();
        TreeModel second = new DefaultTreeModel("second",
                new NewickImporter("((A:0.2,B:0.2):0.1,(C:0.15,(D:0.05,E:0.05):0.1):0.15);").importTree(null));
        List<TreeModel> trees = new ArrayList<>(Arrays.asList(first, second));
        CTMCScalePrior prior = new CTMCScalePrior("prior", new Parameter.Default(0.8), trees, null, false, null,
                false, 0.6);
        check(first, prior, "first of two trees");
        check(second, prior, "second of two trees");
    }

    public void testRefusesOtherTree() throws Exception {
        TreeModel tree = tree();
        CTMCScalePrior prior = new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, false, null, false);
        try {
            prior.getGradientLogDensityWrtNodeHeights(tree());
            fail("a tree the prior does not use must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    public void testReportPasses() throws Exception {
        TreeModel tree = tree();
        CTMCScalePrior prior = new CTMCScalePrior("prior", new Parameter.Default(0.8), tree, null, false, null, false);
        new CTMCScalePriorNodeHeightGradient(prior, tree, 1E-4).getReport(); // throws if the numeric check fails
    }

    private static void check(TreeModel tree, CTMCScalePrior prior, String label) {
        double[] analytic = prior.getGradientLogDensityWrtNodeHeights(tree);
        final double h = 1E-6;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            final double height = tree.getNodeHeight(node);
            tree.setNodeHeight(node, height + h);
            final double plus = prior.getLogLikelihood();
            tree.setNodeHeight(node, height - h);
            final double minus = prior.getLogLikelihood();
            tree.setNodeHeight(node, height);
            final double numeric = (plus - minus) / (2 * h);
            assertEquals(label + " node " + node.getNumber(), numeric, analytic[i],
                    1E-6 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static HKY hky() {
        return new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.1, 0.2, 0.3, 0.4}));
    }

    private static TreeModel tree() throws Exception {
        Tree newick = new NewickImporter(NEWICK).importTree(null);
        return new DefaultTreeModel("treeModel", newick);
    }
}
