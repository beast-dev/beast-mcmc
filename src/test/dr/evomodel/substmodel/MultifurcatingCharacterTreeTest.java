/*
 * MultifurcatingCharacterTreeTest.java
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

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.AbstractHierarchicalRateProvider;
import dr.evomodel.substmodel.CharacterTreeContraction;
import dr.evomodel.substmodel.EqualRatesProvider;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.RelativeRatesProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Multifurcating character trees held as a binary skeleton with a collapse state.
 *
 * @author Xiang Ji
 */
public class MultifurcatingCharacterTreeTest extends TestCase {

    private static final double[] PI = {0.1, 0.2, 0.3, 0.4};

    private static FrequencyModel frequencies(DataType type, double[] pi) {
        return new FrequencyModel(type, new Parameter.Default(pi.clone()));
    }

    private static double maxDiff(SubstitutionModel a, SubstitutionModel b) {
        final int n = a.getDataType().getStateCount();
        double[] x = new double[n * n], y = new double[n * n];
        a.getInfinitesimalMatrix(x);
        b.getInfinitesimalMatrix(y);
        double d = 0.0;
        for (int i = 0; i < n * n; i++) d = Math.max(d, Math.abs(x[i] - y[i]));
        for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
            a.getTransitionProbabilities(t, x);
            b.getTransitionProbabilities(t, y);
            for (int i = 0; i < n * n; i++) d = Math.max(d, Math.abs(x[i] - y[i]));
        }
        return d;
    }

    private static CharacterTreeContraction fixed(String newick, DataType type) {
        CharacterTreeContraction.Built built = CharacterTreeContraction.build("characterTree", newick, type, true, null);
        return new CharacterTreeContraction(built.skeleton, type, built.collapsed, null);
    }

    private static String clade(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return tree.getTaxonId(node.getNumber());
        }
        List<String> codes = new ArrayList<String>();
        for (int i = 0; i < tree.getChildCount(node); i++) {
            codes.add(clade(tree, tree.getChild(node, i)));
        }
        char[] chars = String.join("", codes).toCharArray();
        Arrays.sort(chars);
        return new String(chars);
    }

    /**
     * Raw rates on the skeleton from rates by clade name; collapsed nodes get a value that must not matter.
     */
    private static Parameter skeletonRates(TreeModel skeleton, Map<String, Double> byClade, double collapsedValue) {
        final int tips = skeleton.getExternalNodeCount();
        Parameter rates = new Parameter.Default(skeleton.getInternalNodeCount(), 0.0);
        for (int i = 0; i < skeleton.getInternalNodeCount(); i++) {
            Double r = byClade.get(clade(skeleton, skeleton.getNode(tips + i)));
            rates.setParameterValue(i, r == null ? collapsedValue : r);
        }
        return rates;
    }

    private static HierarchicalSubstitutionModel defaultModel(CharacterTreeContraction contraction, Map<String, Double> byClade,
                                                              FrequencyModel f, double collapsedValue) {
        return new HierarchicalSubstitutionModel("hrate", new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(
                skeletonRates((TreeModel) contraction.getTree(), byClade, collapsedValue), contraction, f, true), f);
    }

    public void testStarIsF81() throws Exception {
        FrequencyModel f = frequencies(Nucleotides.INSTANCE, PI);
        HKY f81 = new HKY(new Parameter.Default(1, 1.0), f);
        CharacterTreeContraction star = fixed("(A,C,G,T);", Nucleotides.INSTANCE);
        assertEquals(1, star.getContractedInternalNodeCount());
        assertEquals("(A,C,G,T);", star.contractedNewick(null));
        assertEquals(Math.log(15.0), star.getLogResolutionCount(), 1e-14);

        HierarchicalSubstitutionModel equal = new HierarchicalSubstitutionModel("hrate", new EqualRatesProvider(star, f, true));
        assertEquals(0.0, maxDiff(equal, f81), 1e-14);
        assertEquals(1, equal.getContractedInternalNodeCount());
        assertEquals("ACGT", String.join(",", ((AbstractHierarchicalRateProvider) equal.getRateProvider()).getCladeNames()));

        Map<String, Double> rootOnly = new HashMap<String, Double>();
        rootOnly.put("ACGT", 2.7);
        HierarchicalSubstitutionModel defaultStar = defaultModel(star, rootOnly, f, 5.0);
        assertEquals(0.0, maxDiff(defaultStar, f81), 1e-14);
        assertTrue(defaultStar.isValidGenerator());

        double[] triples = equal.getBeagleModelTriples();
        double[] lambdas = equal.getBeagleModelLambdas();
        assertEquals("[4.0, 0.0, 1.0, 5.0, 2.0, 3.0, 6.0, 4.0, 5.0]", Arrays.toString(Arrays.copyOf(triples, 9)));
        assertEquals(0.0, lambdas[1], 0.0);
        assertEquals(0.0, lambdas[2], 0.0);
        assertTrue(lambdas[3] > 0.0);
        assertTrue(equal.hasNonNegativeRates());
        double[] other = defaultStar.getBeagleModelLambdas();
        for (int k = 0; k < lambdas.length; k++) assertEquals(lambdas[k], other[k], 1e-15 * Math.abs(lambdas[k]));
    }

    public void testOneCladeTreeIsTN93WithKappa2One() throws Exception {
        FrequencyModel f = frequencies(Nucleotides.INSTANCE, PI);
        final double piAG = PI[0] + PI[2];
        final double rootRate = 1.1, agRate = 0.8;
        TN93 tn93 = new TN93(new Parameter.Default(1, (rootRate + agRate / piAG) / rootRate), new Parameter.Default(1, 1.0), f);
        for (String newick : new String[]{"((A,G),C,T);", "(T,(G,A),C);", "(C,T,(A:2,G:3):1);"}) {
            CharacterTreeContraction contraction = fixed(newick, Nucleotides.INSTANCE);
            assertEquals(2, contraction.getContractedInternalNodeCount());
            assertEquals("((A,G),C,T);", contraction.contractedNewick(null));
            assertEquals(Math.log(3.0), contraction.getLogResolutionCount(), 1e-14);
            Map<String, Double> rates = new HashMap<String, Double>();
            rates.put("ACGT", rootRate);
            rates.put("AG", agRate);
            HierarchicalSubstitutionModel model = defaultModel(contraction, rates, f, 7.0);
            assertEquals(newick, 0.0, maxDiff(model, tn93), 1e-14);
            assertEquals(2, model.getContractedInternalNodeCount());
        }

        CharacterTreeContraction contraction = fixed("((A,G),C,T);", Nucleotides.INSTANCE);
        Parameter tied = new Parameter.Default(new double[]{1.0, 3.0, 1.0, 1.0, 1.0, 1.0});
        HierarchicalSubstitutionModel gtr = new HierarchicalSubstitutionModel("hrate",
                RelativeRatesProvider.gtr(tied, f, contraction, 1e-10), f);
        TN93 kappa3 = new TN93(new Parameter.Default(1, 3.0), new Parameter.Default(1, 1.0), f);
        assertEquals(0.0, maxDiff(gtr, kappa3), 1e-14);
        assertEquals("ACGT,AG", String.join(",", ((AbstractHierarchicalRateProvider) gtr.getRateProvider()).getCladeNames()));

        Parameter untied = new Parameter.Default(new double[]{1.0, 3.0, 1.0, 1.0, 2.0, 1.0});
        try {
            RelativeRatesProvider.gtr(untied, f, fixed("((A,G),C,T);", Nucleotides.INSTANCE), 1e-10);
            fail("CT differs from the root class inside the multifurcation");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("not tree-structured on the multifurcating character tree"));
        }
    }

    private static List<String> allFourStateTrees() {
        String[] s = {"A", "C", "G", "T"};
        Set<String> trees = new HashSet<String>();
        trees.add("(A,C,G,T);");
        for (int a = 0; a < 4; a++) {
            for (int b = a + 1; b < 4; b++) {
                List<String> rest = new ArrayList<String>();
                for (int c = 0; c < 4; c++) {
                    if (c != a && c != b) rest.add(s[c]);
                }
                trees.add("((" + s[a] + "," + s[b] + ")," + rest.get(0) + "," + rest.get(1) + ");");
                trees.add("((" + s[a] + "," + s[b] + "),(" + rest.get(0) + "," + rest.get(1) + "));");
                trees.add("(((" + s[a] + "," + s[b] + ")," + rest.get(0) + ")," + rest.get(1) + ");");
                trees.add("(((" + s[a] + "," + s[b] + ")," + rest.get(1) + ")," + rest.get(0) + ");");
            }
        }
        for (int w = 0; w < 4; w++) {
            List<String> rest = new ArrayList<String>();
            for (int c = 0; c < 4; c++) {
                if (c != w) rest.add(s[c]);
            }
            trees.add("((" + String.join(",", rest) + ")," + s[w] + ");");
        }
        List<String> canonical = new ArrayList<String>();
        Set<String> seen = new HashSet<String>();
        for (String t : trees) {
            String c = fixed(t, Nucleotides.INSTANCE).contractedNewick(null);
            if (seen.add(c)) {
                canonical.add(t);
            }
        }
        return canonical;
    }

    /**
     * Q from the definition on the multifurcating tree as written: Q_ij = pi_j sigma(mrca(i, j)), sigma(m) the sum of
     * rate / pi over the path from m to the root, normalised to one expected substitution per unit time.
     */
    private static double[] definitionQ(String newick, Map<String, Double> rates, double[] pi) throws Exception {
        Tree tree = new NewickImporter(newick).importTree(null);
        final int n = pi.length;
        Map<String, Double> sigma = new HashMap<String, Double>();
        fillSigma(tree, tree.getRoot(), 0.0, rates, pi, sigma);
        double[] q = new double[n * n];
        String codes = "ACGT";
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                String mrca = mrca(tree, tree.getRoot(), codes.charAt(i), codes.charAt(j));
                q[i * n + j] = pi[j] * sigma.get(mrca);
                sum += q[i * n + j];
            }
            q[i * n + i] = -sum;
        }
        double mu = 0.0;
        for (int i = 0; i < n; i++) mu -= pi[i] * q[i * n + i];
        for (int k = 0; k < n * n; k++) q[k] /= mu;
        return q;
    }

    private static void fillSigma(Tree tree, NodeRef node, double above, Map<String, Double> rates, double[] pi,
                                  Map<String, Double> sigma) {
        if (tree.isExternal(node)) return;
        String c = clade(tree, node);
        double piClade = 0.0;
        for (char code : c.toCharArray()) piClade += pi["ACGT".indexOf(code)];
        double s = above + rates.get(c) / piClade;
        sigma.put(c, s);
        for (int i = 0; i < tree.getChildCount(node); i++) fillSigma(tree, tree.getChild(node, i), s, rates, pi, sigma);
    }

    private static String mrca(Tree tree, NodeRef node, char a, char b) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            String c = clade(tree, tree.getChild(node, i));
            if (c.indexOf(a) >= 0 && c.indexOf(b) >= 0) return mrca(tree, tree.getChild(node, i), a, b);
        }
        return clade(tree, node);
    }

    private static double[] expm(double[] q, double t, int n) {
        int squarings = 0;
        double norm = 0.0;
        for (double v : q) norm = Math.max(norm, Math.abs(v * t));
        while (norm > 0.01) {
            norm /= 2.0;
            squarings++;
        }
        final double scale = t / Math.pow(2.0, squarings);
        double[] result = new double[n * n], term = new double[n * n];
        for (int i = 0; i < n; i++) {
            result[i * n + i] = 1.0;
            term[i * n + i] = 1.0;
        }
        for (int k = 1; k <= 18; k++) {
            term = multiply(term, q, n);
            for (int i = 0; i < n * n; i++) term[i] *= scale / k;
            for (int i = 0; i < n * n; i++) result[i] += term[i];
        }
        for (int s = 0; s < squarings; s++) result = multiply(result, result, n);
        return result;
    }

    private static double[] multiply(double[] a, double[] b, int n) {
        double[] c = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                for (int j = 0; j < n; j++) c[i * n + j] += a[i * n + k] * b[k * n + j];
        return c;
    }

    private static void collectClades(Tree tree, NodeRef node, List<String> out) {
        if (tree.isExternal(node)) return;
        out.add(clade(tree, node));
        for (int i = 0; i < tree.getChildCount(node); i++) collectClades(tree, tree.getChild(node, i), out);
    }

    public void testAllFourStateTreesAgainstTheDefinition() throws Exception {
        List<String> trees = allFourStateTrees();
        assertEquals(26, trees.size());
        int[] byQ = new int[3];
        Random random = new Random(20261114);
        FrequencyModel f = frequencies(Nucleotides.INSTANCE, PI);
        for (String newick : trees) {
            CharacterTreeContraction contraction = fixed(newick, Nucleotides.INSTANCE);
            final int q = contraction.getContractedInternalNodeCount() - 1;
            byQ[q]++;
            assertEquals(newick, q == 0 ? Math.log(15.0) : q == 1 ? Math.log(3.0) : 0.0,
                    contraction.getLogResolutionCount(), 1e-14);
            Tree written = new NewickImporter(newick).importTree(null);
            List<String> clades = new ArrayList<String>();
            collectClades(written, written.getRoot(), clades);
            Map<String, Double> rates = new HashMap<String, Double>();
            for (String c : clades) rates.put(c, 0.2 + 2.0 * random.nextDouble());
            HierarchicalSubstitutionModel model = defaultModel(contraction, rates, f, 3.0);
            double[] expected = definitionQ(newick, rates, PI);
            double[] actual = new double[16];
            model.getInfinitesimalMatrix(actual);
            for (int k = 0; k < 16; k++) assertEquals(newick + " Q " + k, expected[k], actual[k], 1e-13);
            for (double t : new double[]{0.01, 0.3, 2.0}) {
                double[] p = expm(expected, t, 4);
                model.getTransitionProbabilities(t, actual);
                for (int k = 0; k < 16; k++) assertEquals(newick + " P(" + t + ") " + k, p[k], actual[k], 1e-12);
                double[] fromEigen = fromEigen(model, t);
                for (int k = 0; k < 16; k++) assertEquals(newick + " eigen P(" + t + ") " + k, p[k], fromEigen[k], 1e-12);
            }
        }
        assertEquals("[1, 10, 15]", Arrays.toString(byQ));
    }

    private static double[] fromEigen(HierarchicalSubstitutionModel model, double t) {
        dr.evomodel.substmodel.EigenDecomposition ed = model.getEigenDecomposition();
        final int n = 4;
        double[] evec = ed.getEigenVectors(), ievc = ed.getInverseEigenVectors(), eval = ed.getEigenValues();
        double[] p = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) sum += evec[i * n + k] * Math.exp(eval[k] * t) * ievc[k * n + j];
                p[i * n + j] = sum;
            }
        return p;
    }

    private static String randomMultifurcatingNewick(List<String> codes, Random random, int maxArity) {
        List<String> pieces = new ArrayList<String>(codes);
        while (pieces.size() > 1) {
            final int k = Math.min(pieces.size(), 2 + random.nextInt(maxArity - 1));
            if (pieces.size() - k == 1 || pieces.size() == k) {
                if (pieces.size() <= maxArity) {
                    return "(" + String.join(",", pieces) + ");";
                }
            }
            List<String> chosen = new ArrayList<String>();
            for (int i = 0; i < k; i++) chosen.add(pieces.remove(random.nextInt(pieces.size())));
            pieces.add(random.nextInt(pieces.size() + 1), "(" + String.join(",", chosen) + ")");
        }
        return pieces.get(0) + ";";
    }

    private static GeneralDataType letters(int n) {
        String[] codes = new String[n];
        for (int i = 0; i < n; i++) codes[i] = String.valueOf((char) ('a' + i));
        return new GeneralDataType(codes);
    }

    public void testUploadDoesNotDependOnTheResolution() throws Exception {
        Random random = new Random(1114);
        int multifurcating = 0;
        for (int trial = 0; trial < 120; trial++) {
            final int n = 4 + random.nextInt(9);
            GeneralDataType type = letters(n);
            List<String> codes = new ArrayList<String>();
            for (int i = 0; i < n; i++) codes.add(type.getCode(i));
            String newick = randomMultifurcatingNewick(codes, random, 5);
            double[] pi = new double[n];
            double total = 0.0;
            for (int i = 0; i < n; i++) total += pi[i] = -Math.log(random.nextDouble());
            for (int i = 0; i < n; i++) pi[i] /= total;
            Tree written = new NewickImporter(newick).importTree(null);
            List<String> clades = new ArrayList<String>();
            collectClades(written, written.getRoot(), clades);
            Map<String, Double> rates = new HashMap<String, Double>();
            for (String c : clades) rates.put(c, 0.1 + random.nextDouble());
            double[] triples = null, lambdas = null;
            String contracted = null;
            for (int r = 0; r < 6; r++) {
                CharacterTreeContraction.Built built = CharacterTreeContraction.build("characterTree", newick, type, true,
                        r == 0 ? null : random);
                multifurcating += r == 0 && built.isMultifurcating() ? 1 : 0;
                CharacterTreeContraction contraction = new CharacterTreeContraction(built.skeleton, type, built.collapsed, null);
                FrequencyModel f = frequencies(type, pi);
                HierarchicalSubstitutionModel model = defaultModel(contraction, rates, f, 0.5 + r);
                if (triples == null) {
                    triples = model.getBeagleModelTriples().clone();
                    lambdas = model.getBeagleModelLambdas().clone();
                    contracted = contraction.contractedNewick(null);
                } else {
                    assertEquals(newick, contracted, contraction.contractedNewick(null));
                    assertTrue(newick + " triples", Arrays.equals(triples, model.getBeagleModelTriples()));
                    assertTrue(newick + " lambdas resolution " + r, Arrays.equals(lambdas, model.getBeagleModelLambdas()));
                }
            }
        }
        assertTrue(multifurcating > 60);
    }

    private static final class Counter implements ModelListener {
        int fired = 0;

        public void modelChangedEvent(Model model, Object object, int index) {
            fired++;
        }

        public void modelRestored(Model model) {
        }
    }

    private static NodeRef nodeOfClade(TreeModel tree, String clade) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            if (clade(tree, tree.getNode(i)).equals(clade)) return tree.getNode(i);
        }
        throw new IllegalArgumentException(clade);
    }

    private static void exchange(TreeModel tree, NodeRef child, NodeRef other) {
        NodeRef p1 = tree.getParent(child), p2 = tree.getParent(other);
        tree.beginTreeEdit();
        tree.removeChild(p1, child);
        tree.removeChild(p2, other);
        tree.addChild(p1, other);
        tree.addChild(p2, child);
        tree.endTreeEdit();
    }

    public void testNeutralMovesFireNothing() throws Exception {
        GeneralDataType type = letters(6);
        FrequencyModel f = frequencies(type, new double[]{0.07, 0.11, 0.13, 0.17, 0.23, 0.29});
        String newick = "((a,b,c,d),(e,f));";
        CharacterTreeContraction.Built built = CharacterTreeContraction.build("characterTree", newick, type, false, null);
        TreeModel skeleton = built.skeleton;
        assertTrue(contraction(skeleton, "ab") && contraction(skeleton, "cd"));
        CharacterTreeContraction contraction = new CharacterTreeContraction(skeleton, type, built.collapsed, null);
        Map<String, Double> rates = new HashMap<String, Double>();
        rates.put("abcdef", 1.0);
        rates.put("abcd", 0.4);
        rates.put("ef", 2.0);
        Parameter raw = skeletonRates(skeleton, rates, 9.0);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(raw, contraction, f, true), f);
        Counter counter = new Counter();
        model.addModelListener(counter);
        final double[] triples = model.getBeagleModelTriples().clone();
        final double[] lambdas = model.getBeagleModelLambdas().clone();
        double[] before = new double[36];
        model.getTransitionProbabilities(0.4, before);

        model.storeModelState();
        exchange(skeleton, nodeOfClade(skeleton, "b"), nodeOfClade(skeleton, "c"));
        assertEquals("a move inside the collapsed cluster abcd", 0, counter.fired);
        assertEquals("((a,b,c,d),(e,f));", contraction.contractedNewick(null));
        assertTrue(Arrays.equals(triples, model.getBeagleModelTriples()));
        assertTrue(Arrays.equals(lambdas, model.getBeagleModelLambdas()));
        double[] after = new double[36];
        model.getTransitionProbabilities(0.4, after);
        for (int k = 0; k < 36; k++) assertEquals(before[k], after[k], 1e-15);

        HierarchicalSubstitutionModel.HierarchicalRateProvider.Default provider =
                (HierarchicalSubstitutionModel.HierarchicalRateProvider.Default) model.getRateProvider();
        NodeRef collapsed = null;
        for (int i = skeleton.getExternalNodeCount(); i < skeleton.getNodeCount(); i++) {
            if (contraction.isCollapsed(i)) collapsed = skeleton.getNode(i);
        }
        assertNotNull(collapsed);
        raw.setParameterValue(collapsed.getNumber() - skeleton.getExternalNodeCount(), 42.0);
        assertEquals("a change of a collapsed rate", 0, counter.fired);
        assertTrue(Arrays.equals(lambdas, model.getBeagleModelLambdas()));
        assertTrue(provider.isNormalised());

        exchange(skeleton, nodeOfClade(skeleton, "a"), nodeOfClade(skeleton, "e"));
        assertTrue("a move that changes the contracted tree", counter.fired > 0);
        assertFalse(Arrays.equals(triples, model.getBeagleModelTriples()) && Arrays.equals(lambdas, model.getBeagleModelLambdas()));

        model.restoreModelState();
        assertEquals(9.0, raw.getParameterValue(collapsed.getNumber() - skeleton.getExternalNodeCount()), 0.0);
        assertEquals("((a,b,c,d),(e,f));", contraction.contractedNewick(null));
        final int firedBefore = counter.fired;
        exchange(skeleton, nodeOfClade(skeleton, "a"), nodeOfClade(skeleton, "d"));
        assertEquals("a neutral move after the restore", firedBefore, counter.fired);
        assertTrue(Arrays.equals(triples, model.getBeagleModelTriples()));
        assertTrue(Arrays.equals(lambdas, model.getBeagleModelLambdas()));
    }

    private static boolean contraction(TreeModel skeleton, String clade) {
        for (int i = skeleton.getExternalNodeCount(); i < skeleton.getNodeCount(); i++) {
            if (clade(skeleton, skeleton.getNode(i)).equals(clade)) return true;
        }
        return false;
    }

    public void testKeepParameterContraction() throws Exception {
        TreeModel skeleton = AbstractHierarchicalRateProvider.characterTree("(((A,C),G),T);");
        Parameter keep = new Parameter.Default(3, 1.0);
        CharacterTreeContraction contraction = new CharacterTreeContraction(skeleton, Nucleotides.INSTANCE, null, keep);
        assertEquals("(((A,C),G),T);", contraction.contractedNewick(null));
        assertEquals(0.0, contraction.getLogResolutionCount(), 0.0);
        int ac = nodeOfClade(skeleton, "AC").getNumber(), acg = nodeOfClade(skeleton, "ACG").getNumber();
        int root = skeleton.getRoot().getNumber();
        keep.setParameterValue(ac - 4, 0.0);
        assertEquals("((A,C,G),T);", contraction.contractedNewick(null));
        assertEquals(Math.log(3.0), contraction.getLogResolutionCount(), 1e-15);
        keep.setParameterValue(acg - 4, 0.0);
        assertEquals("(A,C,G,T);", contraction.contractedNewick(null));
        assertEquals(1, contraction.getContractedInternalNodeCount());
        keep.setParameterValue(root - 4, 0.0);
        assertFalse("the root is never collapsed", contraction.isCollapsed(root));
        assertEquals("(A,C,G,T);", contraction.contractedNewick(null));

        FrequencyModel f = frequencies(Nucleotides.INSTANCE, PI);
        Parameter rates = new Parameter.Default(new double[]{0.3, 0.6, 1.0});
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, contraction, f, true), f);
        assertEquals(0.0, maxDiff(model, new HKY(new Parameter.Default(1, 1.0), f)), 1e-14);
        assertEquals(0.0, rates.getBounds().getLowerLimit(0), 0.0);
        Counter counter = new Counter();
        model.addModelListener(counter);
        keep.setParameterValue(acg - 4, 1.0);
        assertTrue("a change of keep fires", counter.fired > 0);
        assertEquals("((A,C,G),T);", contraction.contractedNewick(null));
        try {
            new EqualRatesProvider(contraction, f, true);
            fail("a keep parameter needs ratesParameter");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("needs ratesParameter"));
        }
    }

    public void testEqualRatesOnRandomMultifurcatingTrees() throws Exception {
        Random random = new Random(7);
        for (int n : new int[]{2, 4, 7, 20}) {
            GeneralDataType type = letters(n);
            List<String> codes = new ArrayList<String>();
            for (int i = 0; i < n; i++) codes.add(type.getCode(i));
            double[] pi = new double[n];
            double total = 0.0;
            for (int i = 0; i < n; i++) total += pi[i] = 0.2 + random.nextDouble();
            double sumSquares = 0.0;
            for (int i = 0; i < n; i++) {
                pi[i] /= total;
                sumSquares += pi[i] * pi[i];
            }
            for (int trial = 0; trial < 5; trial++) {
                String newick = n == 2 ? "(a,b);" : trial == 0 ? CharacterTreeContraction.starNewick(type)
                        : randomMultifurcatingNewick(codes, random, 4);
                CharacterTreeContraction.Built built = CharacterTreeContraction.build("characterTree", newick, type, true, null);
                FrequencyModel f = frequencies(type, pi);
                HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", built.isMultifurcating()
                        ? new EqualRatesProvider(new CharacterTreeContraction(built.skeleton, type, built.collapsed, null), f, true)
                        : new EqualRatesProvider(built.skeleton, f, true));
                double[] q = new double[n * n];
                model.getInfinitesimalMatrix(q);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            assertEquals(newick, pi[j] / (1.0 - sumSquares), q[i * n + j], 1e-13);
                        }
                    }
                }
                double[] p = new double[n * n];
                model.getTransitionProbabilities(0.7, p);
                final double e = Math.exp(-0.7 / (1.0 - sumSquares));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        assertEquals(newick, (i == j ? e : 0.0) + (1.0 - e) * pi[j], p[i * n + j], 1e-13);
                    }
                }
            }
        }
    }
}
