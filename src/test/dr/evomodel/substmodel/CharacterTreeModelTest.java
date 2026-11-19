/*
 * CharacterTreeModelTest.java
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

package test.dr.evomodel.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evomodel.substmodel.CharacterTreeStatistic;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.substmodel.CharacterTreeModelParser;
import dr.evomodelxml.substmodel.CharacterTreeStatisticParser;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.HierarchicalSubstitutionModelParser;
import dr.evoxml.NewickParser;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The character-tree model: its contracted view, events, counts, guard and parsing.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeModelTest extends TestCase {

    private static final class Counter implements ModelListener {
        int fired = 0;

        public void modelChangedEvent(Model model, Object object, int index) {
            fired++;
        }

        public void modelRestored(Model model) {
        }
    }

    private static NodeRef node(DefaultTreeModel tree, String clade) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            if (clade(tree, tree.getNode(i)).equals(clade)) return tree.getNode(i);
        }
        throw new IllegalArgumentException(clade);
    }

    private static String clade(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) return tree.getTaxonId(node.getNumber());
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getChildCount(node); i++) sb.append(clade(tree, tree.getChild(node, i)));
        char[] c = sb.toString().toCharArray();
        java.util.Arrays.sort(c);
        return new String(c);
    }

    private static void exchangeWithUncle(CharacterTreeModel model, NodeRef u, int which) {
        DefaultTreeModel tree = model.getSkeleton();
        NodeRef p = tree.getParent(u);
        NodeRef w = tree.getChild(p, 0) == u ? tree.getChild(p, 1) : tree.getChild(p, 0);
        NodeRef c = tree.getChild(u, which);
        tree.beginTreeEdit();
        tree.removeChild(u, c);
        tree.removeChild(p, w);
        tree.addChild(u, w);
        tree.addChild(p, c);
        model.resetUnitHeightsQuietly();
        tree.endTreeEdit();
    }

    public void testViewIsTheContractedTree() throws Exception {
        Parameter keep = new Parameter.Default(1);
        CharacterTreeModel model = new CharacterTreeModel("characterTree", "((T,G),C,A);", Nucleotides.INSTANCE, true, keep);
        assertEquals(3, keep.getDimension());
        assertEquals("(A,C,(G,T));", model.getContractedNewick());
        assertEquals(2, model.getContractedInternalNodeCount());
        assertEquals(Math.log(3.0), model.getLogResolutionCount(), 1e-15);
        assertEquals(4, model.getExternalNodeCount());
        assertEquals(2, model.getInternalNodeCount());
        assertEquals(6, model.getNodeCount());
        for (int i = 0; i < 4; i++) {
            NodeRef tip = model.getExternalNode(i);
            assertEquals(i, tip.getNumber());
            assertEquals(model.getSkeleton().getTaxonId(i), model.getNodeTaxon(tip).getId());
        }
        NodeRef root = model.getRoot();
        assertEquals(5, root.getNumber());
        assertEquals(3, model.getChildCount(root));
        assertEquals("A", model.getNodeTaxon(model.getChild(root, 0)).getId());
        assertEquals("C", model.getNodeTaxon(model.getChild(root, 1)).getId());
        assertEquals(4, model.getChild(root, 2).getNumber());
        for (int i = 0; i < model.getNodeCount(); i++) {
            if (!model.isRoot(model.getNode(i))) assertEquals(1.0, model.getBranchLength(model.getNode(i)), 0.0);
        }
        assertEquals("(A:1.0,C:1.0,(G:1.0,T:1.0):1.0);", TreeUtils.newick(model));
        assertEquals(model.getSkeleton().getTaxonCount(), model.getTaxonCount());
    }

    public void testEventsOnlyWhenTheContractedTreeChanges() throws Exception {
        Parameter keep = new Parameter.Default(1);
        CharacterTreeModel model = new CharacterTreeModel("characterTree", "((A,C,G),T);", Nucleotides.INSTANCE, true, keep);
        Counter counter = new Counter();
        model.addModelListener(counter);
        DefaultTreeModel skeleton = model.getSkeleton();
        NodeRef inserted = null;
        for (int i = 4; i < 7; i++) {
            if (model.isCollapsed(i)) inserted = skeleton.getNode(i);
        }
        assertNotNull(inserted);
        exchangeWithUncle(model, inserted, 0);
        assertEquals("a move inside the collapsed cluster", 0, counter.fired);
        assertEquals("((A,C,G),T);", model.getContractedNewick());

        keep.setParameterValue(inserted.getNumber() - 4, 1.0);
        assertTrue(counter.fired > 0);
        assertFalse("((A,C,G),T);".equals(model.getContractedNewick()));
        final int before = counter.fired;
        NodeRef acg = node(skeleton, "ACG");
        exchangeWithUncle(model, acg, 0);
        assertTrue("a move of a kept node", counter.fired > before);
    }

    public void testCountsOfTheFourStateTrees() throws Exception {
        Map<String, Double> logR = new HashMap<String, Double>();
        int[] byQ = new int[3];
        List<String> trees = CharacterTreeStatistic.treeList(Nucleotides.INSTANCE);
        assertEquals(26, trees.size());
        for (String t : trees) {
            CharacterTreeModel model = new CharacterTreeModel("characterTree", t, Nucleotides.INSTANCE, true, new Parameter.Default(1));
            assertEquals(t, model.getContractedNewick());
            final int q = model.getContractedInternalNodeCount() - 1;
            byQ[q]++;
            logR.put(t, model.getLogResolutionCount());
            assertEquals(t, q == 0 ? Math.log(15) : q == 1 ? Math.log(3) : 0.0, model.getLogResolutionCount(), 1e-14);
        }
        assertEquals("[1, 10, 15]", java.util.Arrays.toString(byQ));
        double sumR = 0.0;
        for (double v : logR.values()) sumR += Math.exp(v);
        assertEquals("15 + 10 x 3 + 15", 60.0, sumR, 1e-12);
    }

    public void testGuardRefusesGenericTreeOperators() throws Exception {
        MathUtils.setSeed(20261114);
        for (int which = 0; which < 2; which++) {
            CharacterTreeModel model = new CharacterTreeModel("characterTree", "(((A,C),G),T);", Nucleotides.INSTANCE, true,
                    new Parameter.Default(1));
            boolean refused = false;
            for (int k = 0; k < 50 && !refused; k++) {
                try {
                    if (which == 0) {
                        new SubtreeSlideOperator(model.getSkeleton(), 1.0, 0.5, false, false, false, false,
                                AdaptationMode.ADAPTATION_OFF, 0.23).doOperation();
                    } else {
                        new WilsonBalding(model.getSkeleton(), 1.0).doOperation();
                    }
                } catch (IllegalStateException e) {
                    assertTrue(e.getMessage(), e.getMessage().contains("can only be changed by the character-tree operators"));
                    refused = true;
                } catch (RuntimeException e) {
                    continue;
                }
            }
            assertTrue(which == 0 ? "subtree slide" : "Wilson-Balding", refused);
        }
    }

    private static Map<String, XMLObject> parse(String body) throws Exception {
        XMLParser parser = new XMLParser(true, true, false, null);
        parser.addXMLObjectParser(new NewickParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new FrequencyModelParser());
        parser.addXMLObjectParser(new CharacterTreeModelParser());
        parser.addXMLObjectParser(new CharacterTreeStatisticParser());
        parser.addXMLObjectParser(new HierarchicalSubstitutionModelParser());
        return parser.parse(new StringReader("<?xml version=\"1.0\" standalone=\"yes\"?>\n<beast>\n" + body + "</beast>\n"), false);
    }

    private static String failure(String body) throws Exception {
        try {
            parse(body);
        } catch (XMLParseException e) {
            return e.getMessage();
        }
        fail("the parse should fail");
        return null;
    }

    private static final String F = "  <frequencyModel id=\"f\" dataType=\"nucleotide\"><frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel>\n";

    private static String hrate(String body) {
        return "  <hierarchicalSubstitutionModel id=\"hrate\">\n    <frequencies><frequencyModel idref=\"f\"/></frequencies>\n"
                + "    <characterTree><characterTreeModel idref=\"characterTree\"/></characterTree>\n" + body
                + "  </hierarchicalSubstitutionModel>\n";
    }

    public void testParsing() throws Exception {
        String search = "  <characterTreeModel id=\"characterTree\" variable=\"true\" newick=\"((A,G),(C,T));\">\n"
                + "    <keep><parameter id=\"characterTree.keep\"/></keep>\n  </characterTreeModel>\n";
        Map<String, XMLObject> store = parse(F + search
                + hrate("    <ratesParameter><parameter id=\"hrate.rates\" value=\"1.0\" lower=\"0.0\"/></ratesParameter>\n"));
        CharacterTreeModel tree = (CharacterTreeModel) store.get("characterTree").getNativeObject();
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        Parameter rates = (Parameter) store.get("hrate.rates").getNativeObject();
        assertEquals(3, rates.getDimension());
        assertEquals(0.0, rates.getBounds().getLowerLimit(0), 0.0);
        assertEquals(3, tree.getKeep().getDimension());
        assertEquals("((A,G),(C,T));", tree.getContractedNewick());
        assertSame(tree.getContraction(), model.getContraction());

        String message = failure(F + search + hrate("    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("must be given with ratesParameter"));
        message = failure(F + search + hrate("    <equalRates/>\n"));
        assertTrue(message, message.contains("must be given with ratesParameter"));
        message = failure(F + "  <characterTreeModel id=\"characterTree\" variable=\"true\" newick=\"((A,G),C,T);\"/>\n");
        assertTrue(message, message.contains("needs keep"));
        message = failure(F + "  <characterTreeModel id=\"characterTree\" variable=\"false\" newick=\"((A,G),C,T);\">"
                + "<keep><parameter value=\"1\"/></keep></characterTreeModel>\n");
        assertTrue(message, message.contains("variable"));
        message = failure("  <characterTreeModel id=\"characterTree\" variable=\"true\" dataType=\"binary\" newick=\"(0,1);\">"
                + "<keep><parameter value=\"1\"/></keep></characterTreeModel>\n");
        assertTrue(message, message.contains("fewer than 3 states"));

        store = parse(F + "  <characterTreeModel id=\"characterTree\" newick=\"((A,G),C,T);\"/>\n"
                + hrate("    <ratesParameter><parameter id=\"hrate.rates\" value=\"0.8 1.1\"/></ratesParameter>\n"));
        model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        assertTrue(model.getRateProvider() instanceof dr.evomodel.substmodel.ContractedRatesProvider);
        message = failure(F + "  <characterTreeModel id=\"characterTree\" newick=\"((A,G),C,T);\"/>\n"
                + hrate("    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("needs the binary character tree"));
    }

    public void testRateTraitOnTheView() throws Exception {
        for (String newick : new String[]{"((A,G),C,T);", "((A,C),(G,T));"}) {
            Map<String, XMLObject> store = parse(F + "  <characterTreeModel id=\"characterTree\" newick=\"" + newick + "\"/>\n"
                    + hrate("    <ratesParameter><parameter id=\"hrate.rates\" value=\""
                    + (newick.startsWith("((A,G),C") ? "0.8 1.1" : "0.6 2.1 0.9") + "\"/></ratesParameter>\n"));
            CharacterTreeModel tree = (CharacterTreeModel) store.get("characterTree").getNativeObject();
            HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
            String logged = TreeUtils.newick(tree, new TreeTraitProvider[]{tree});
            Matcher m = Pattern.compile("\\)\\[&rate=([^\\]]+)\\]").matcher(logged);
            int count = 0;
            while (m.find()) count++;
            assertEquals(logged, tree.getInternalNodeCount(), count);
            for (int i = tree.getExternalNodeCount(); i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                double rate = (Double) tree.getTreeTrait("rate").getTrait(tree, node);
                String clade = clade(tree, node);
                NodeRef kept = node(tree.getSkeleton(), clade);
                assertEquals(clade, model.getCladeRate(kept.getNumber()), rate, 0.0);
                assertTrue(rate > 0.0);
            }
            assertTrue(logged, logged.contains("[&rate="));
        }
    }
}
