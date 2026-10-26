/*
 * NodeHeightGradientAnyClockTest.java
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

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightGradientForDiscreteTrait;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * Node-height gradients of the tree likelihood (plan 3, B8): for a clock that is not a DifferentiableBranchRates (they
 * only need the branch rates), summed over the partitions of a compound likelihood, and a clear error without
 * usePreOrder; checked against finite differences.
 */
public class NodeHeightGradientAnyClockTest extends TestCase {

    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final String[][] SEQUENCES = {
            {"A", "ACGTACGTTACGATCGATCGTAGCTAGCTAGCATCG"},
            {"B", "ACGTACGTAACGATCGTTCGTAGCTAGCTAGCATCA"},
            {"C", "ACGAACGTTACGTTCGATCGTAGGTAGCTAGGATCG"},
            {"D", "TCGTACGATACGATCGATCCTAGCTAGCAAGCATCG"},
            {"E", "TCGTACGATACGATCGATCCTAGCTTGCAAGCATGG"},
    };

    public NodeHeightGradientAnyClockTest(String name) {
        super(name);
    }

    // a fixed rate per branch, with no derivatives
    private static class FixedBranchRates extends AbstractBranchRateModel {
        FixedBranchRates() {
            super("fixedBranchRates");
        }

        public double getBranchRate(Tree tree, NodeRef node) {
            return 0.6 + 0.15 * node.getNumber();
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) { }

        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { }

        protected void storeState() { }

        protected void restoreState() { }

        protected void acceptState() { }
    }

    public void testClockWithoutRateDerivatives() throws Exception {
        TreeModel tree = new DefaultTreeModel("treeModel", new NewickImporter(NEWICK).importTree(null));
        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        for (String[] row : SEQUENCES) {
            alignment.addSequence(new Sequence(tree.getTaxon(tree.getTaxonIndex(row[0])), row[1]));
        }
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                new SitePatterns(alignment, null, 0, -1, 1, true),
                new HomogeneousBranchModel(new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.1, 0.2, 0.3, 0.4}))),
                new GammaSiteRateModel("siteModel"), true, false, PartialsRescalingScheme.ALWAYS, false,
                new PreOrderSettings(true, false, true, true));
        TreeDataLikelihood likelihood = new TreeDataLikelihood(delegate, tree, new FixedBranchRates());

        check(new NodeHeightGradientForDiscreteTrait("Sequence", likelihood, delegate, null), likelihood, tree,
                "fixed branch rates");
    }

    public void testTwoPartitionsOnOneTree() throws Exception {
        Map<String, XMLObject> store = parse(xml(true));
        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider) store.get("gradient").getNativeObject();
        Likelihood joint = (Likelihood) store.get("likelihood").getNativeObject();
        check(gradient, joint, (TreeModel) store.get("treeModel").getNativeObject(), "two partitions");
    }

    public void testUsePreOrderIsRequired() {
        try {
            parse(xml(false));
            fail("a tree likelihood without usePreOrder must be refused");
        } catch (Exception e) {
            String message = e.toString();
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                message += " / " + cause;
            }
            assertTrue(message, message.contains("usePreOrder=\"true\""));
        }
    }

    private static void check(GradientWrtParameterProvider gradient, Likelihood likelihood, TreeModel tree,
                              String label) {
        double[] analytic = gradient.getGradientLogDensity();
        assertEquals(tree.getInternalNodeCount(), analytic.length);
        final double h = 1E-6;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            final double height = tree.getNodeHeight(node);
            tree.setNodeHeight(node, height + h);
            final double plus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height - h);
            final double minus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height);
            final double numeric = (plus - minus) / (2 * h);
            assertEquals(label + " node " + node.getNumber(), numeric, analytic[i],
                    1E-5 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static Map<String, XMLObject> parse(String xml) throws Exception {
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    private static String xml(boolean usePreOrder) {
        StringBuilder sb = new StringBuilder("<beast>\n  <taxa id=\"taxa\">");
        for (String[] row : SEQUENCES) {
            sb.append("<taxon id=\"").append(row[0]).append("\"/>");
        }
        sb.append("</taxa>\n  <alignment id=\"alignment\" dataType=\"nucleotide\">\n");
        for (String[] row : SEQUENCES) {
            sb.append("    <sequence><taxon idref=\"").append(row[0]).append("\"/>").append(row[1]).append("</sequence>\n");
        }
        sb.append("  </alignment>\n");
        sb.append("  <patterns id=\"first\" from=\"1\" to=\"18\"><alignment idref=\"alignment\"/></patterns>\n");
        sb.append("  <patterns id=\"second\" from=\"19\"><alignment idref=\"alignment\"/></patterns>\n");
        sb.append("  <newick id=\"startingTree\">").append(NEWICK).append("</newick>\n");
        sb.append("  <treeModel id=\"treeModel\"><tree idref=\"startingTree\"/>\n" +
                "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
                "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
                "  </treeModel>\n");
        sb.append("  <strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"0.8\"/></rate>" +
                "</strictClockBranchRates>\n");
        for (String id : new String[]{"first", "second"}) {
            String kappa = id.equals("first") ? "2.0" : "5.0";
            sb.append("  <HKYModel id=\"hky.").append(id).append("\"><frequencies><frequencyModel dataType=\"nucleotide\">" +
                    "<frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel></frequencies>" +
                    "<kappa><parameter value=\"").append(kappa).append("\"/></kappa></HKYModel>\n");
            sb.append("  <siteModel id=\"siteModel.").append(id).append("\"><substitutionModel><HKYModel idref=\"hky.")
                    .append(id).append("\"/></substitutionModel></siteModel>\n");
            sb.append("  <treeDataLikelihood id=\"treeLikelihood.").append(id).append("\" usePreOrder=\"")
                    .append(usePreOrder).append("\"><partition><patterns idref=\"").append(id)
                    .append("\"/><siteModel idref=\"siteModel.").append(id).append("\"/></partition>" +
                            "<treeModel idref=\"treeModel\"/><strictClockBranchRates idref=\"branchRates\"/>" +
                            "</treeDataLikelihood>\n");
        }
        sb.append("  <likelihood id=\"likelihood\"><treeDataLikelihood idref=\"treeLikelihood.first\"/>" +
                "<treeDataLikelihood idref=\"treeLikelihood.second\"/></likelihood>\n");
        sb.append("  <nodeHeightGradient id=\"gradient\" traitName=\"Sequence\"><likelihood idref=\"likelihood\"/>" +
                "</nodeHeightGradient>\n");
        sb.append("</beast>\n");
        return sb.toString();
    }
}
