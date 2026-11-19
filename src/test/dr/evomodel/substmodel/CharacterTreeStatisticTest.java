/*
 * CharacterTreeStatisticTest.java
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
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evomodel.substmodel.CharacterTreeStatistic;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * The columns of characterTreeStatistic on the 26 four-state trees.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeStatisticTest extends TestCase {

    private static Set<String> clades(String newick) {
        Set<String> out = new HashSet<String>();
        List<StringBuilder> open = new ArrayList<StringBuilder>();
        String body = newick.substring(0, newick.length() - 1);
        for (int i = 0; i < body.length(); i++) {
            char c = body.charAt(i);
            if (c == '(') {
                open.add(new StringBuilder());
            } else if (c == ')') {
                StringBuilder done = open.remove(open.size() - 1);
                char[] chars = done.toString().toCharArray();
                java.util.Arrays.sort(chars);
                if (!open.isEmpty()) {
                    out.add(new String(chars));
                }
            } else if (Character.isLetter(c)) {
                for (StringBuilder sb : open) sb.append(c);
            }
        }
        return out;
    }

    public void testColumnsOnEveryFourStateTree() throws Exception {
        List<String> trees = CharacterTreeStatistic.treeList(Nucleotides.INSTANCE);
        assertEquals(26, trees.size());
        assertEquals("(A,C,G,T);", trees.get(0));
        Random random = new Random(42);
        double[] pi = {0.1, 0.2, 0.3, 0.4};
        for (int index = 0; index < trees.size(); index++) {
            String newick = trees.get(index);
            Parameter keep = new Parameter.Default(1);
            CharacterTreeModel tree = new CharacterTreeModel("characterTree", newick, Nucleotides.INSTANCE, true, keep);
            double[] raw = new double[3];
            for (int i = 0; i < 3; i++) raw[i] = 0.2 + random.nextDouble();
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(pi.clone()));
            HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                    new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(raw),
                            tree.getContraction(), f, true), f);
            tree.setRateModel(model);
            CharacterTreeStatistic statistic = new CharacterTreeStatistic("ct", tree, model, true, true);
            assertEquals(2 + 1 + 10 + 6, statistic.getDimension());
            assertEquals("ct.q", statistic.getDimensionName(0));
            assertEquals("ct.isStar", statistic.getDimensionName(1));
            assertEquals("ct.model", statistic.getDimensionName(2));
            assertEquals("ct.clade.AC", statistic.getDimensionName(3));
            assertEquals("ct.clade.GT", statistic.getDimensionName(8));
            assertEquals("ct.clade.ACG", statistic.getDimensionName(9));
            assertEquals("ct.clade.CGT", statistic.getDimensionName(12));
            assertEquals("ct.sigma.AC", statistic.getDimensionName(13));
            assertEquals("ct.sigma.GT", statistic.getDimensionName(18));
            Set<String> expected = clades(newick);
            assertEquals(newick, expected.size(), statistic.getStatisticValue(0), 0.0);
            assertEquals(newick, expected.isEmpty() ? 1.0 : 0.0, statistic.getStatisticValue(1), 0.0);
            assertEquals(newick, index, statistic.getStatisticValue(2), 0.0);
            for (int k = 3; k < 13; k++) {
                String clade = statistic.getDimensionName(k).substring("ct.clade.".length());
                assertEquals(newick + " " + clade, expected.contains(clade) ? 1.0 : 0.0, statistic.getStatisticValue(k), 0.0);
            }
            double[] q = new double[16];
            model.getInfinitesimalMatrix(q);
            int k = 13;
            for (int i = 0; i < 4; i++) {
                for (int j = i + 1; j < 4; j++) {
                    assertEquals(newick, q[i * 4 + j] / pi[j], statistic.getStatisticValue(k++), 1e-12);
                }
            }
        }
    }
}
