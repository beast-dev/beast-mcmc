/*
 * CharacterTreeStatistic.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Columns that summarise a sampled character tree: q (the number of non-root internal nodes of the contracted tree),
 * isStar (1 when q = 0: F81, JC69 or Mk), at S <= 5 the index of the tree in the canonical list of all
 * multifurcating trees (getTreeList), at S <= 8 one indicator per possible clade (named by its sorted codes) and, by
 * default at S <= 20, the normalised exchangeabilities sigma(mrca(i, j)) = Q_ij / pi_j, defined on every tree.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeStatistic extends Statistic.Abstract {

    public static final int MAX_STATES_FOR_INDEX = 5;
    public static final int MAX_STATES_FOR_CLADES = 8;
    public static final int MAX_STATES_FOR_EXCHANGEABILITIES = 20;

    private final CharacterTreeModel tree;
    private final HierarchicalSubstitutionModel model;
    private final int stateCount;
    private final List<String> treeList;
    private final Map<String, Integer> indexOfTree;
    private final long[] cladeMasks;
    private final Map<Long, Integer> columnOfClade;
    private final boolean exchangeabilities;
    private final List<String> names = new ArrayList<String>();
    private final double[] values;
    private final double[] q;
    private boolean known = false;

    public CharacterTreeStatistic(String name, CharacterTreeModel tree, HierarchicalSubstitutionModel model,
                                  boolean clades, boolean exchangeabilities) {
        super(name);
        this.tree = tree;
        this.model = model;
        this.stateCount = tree.getExternalNodeCount();
        if (model != null && model.getRateProvider().getTreeModel() != tree.getSkeleton()) {
            throw new IllegalArgumentException("the hierarchical substitution model must be the one of the character tree "
                    + tree.getId());
        }
        if (clades && stateCount > MAX_STATES_FOR_CLADES) {
            throw new IllegalArgumentException("clade indicators are logged for at most " + MAX_STATES_FOR_CLADES + " states");
        }
        if (exchangeabilities && model == null) {
            throw new IllegalArgumentException("the exchangeabilities need the hierarchical substitution model");
        }
        final String prefix = name == null ? "characterTree" : name;
        names.add(prefix + ".q");
        names.add(prefix + ".isStar");
        if (stateCount <= MAX_STATES_FOR_INDEX) {
            treeList = treeList(tree.getDataType());
            indexOfTree = new HashMap<String, Integer>();
            for (int i = 0; i < treeList.size(); i++) {
                indexOfTree.put(treeList.get(i), i);
            }
            names.add(prefix + ".model");
        } else {
            treeList = null;
            indexOfTree = null;
        }
        if (clades) {
            final List<Long> masks = new ArrayList<Long>();
            for (long mask = 1; mask < (1L << stateCount) - 1; mask++) {
                if (Long.bitCount(mask) >= 2) {
                    masks.add(mask);
                }
            }
            Collections.sort(masks, new Comparator<Long>() {
                public int compare(Long a, Long b) {
                    final int sizes = Integer.compare(Long.bitCount(a), Long.bitCount(b));
                    return sizes != 0 ? sizes : codes(a).compareTo(codes(b));
                }
            });
            cladeMasks = new long[masks.size()];
            columnOfClade = new HashMap<Long, Integer>();
            for (int i = 0; i < cladeMasks.length; i++) {
                cladeMasks[i] = masks.get(i);
                columnOfClade.put(cladeMasks[i], names.size());
                names.add(prefix + ".clade." + codes(cladeMasks[i]));
            }
        } else {
            cladeMasks = new long[0];
            columnOfClade = new HashMap<Long, Integer>();
        }
        this.exchangeabilities = exchangeabilities;
        if (exchangeabilities) {
            for (int i = 0; i < stateCount; i++) {
                for (int j = i + 1; j < stateCount; j++) {
                    names.add(prefix + ".sigma." + code(i) + code(j));
                }
            }
        }
        this.values = new double[names.size()];
        this.q = exchangeabilities ? new double[stateCount * stateCount] : null;
    }

    private String code(int state) {
        return AbstractHierarchicalRateProvider.codeOfState(tree.getDataType(), state);
    }

    private String codes(long mask) {
        final List<String> codes = new ArrayList<String>();
        for (int s = 0; s < stateCount; s++) {
            if ((mask & (1L << s)) != 0) {
                codes.add(code(s));
            }
        }
        Collections.sort(codes);
        final StringBuilder sb = new StringBuilder();
        for (String c : codes) {
            sb.append(c);
        }
        return sb.toString();
    }

    /**
     * Every multifurcating tree on the states of a data type (S <= 5), as the contracted Newick strings of
     * CharacterTreeModel, ordered by q and then by the string: the meaning of the model column.
     */
    public static List<String> treeList(DataType dataType) {
        final int n = dataType.getStateCount();
        if (n > MAX_STATES_FOR_INDEX) {
            throw new IllegalArgumentException("the tree list is built for at most " + MAX_STATES_FOR_INDEX + " states");
        }
        final List<String> labels = new ArrayList<String>();
        for (int i = 0; i < n; i++) {
            final String code = AbstractHierarchicalRateProvider.codeOfState(dataType, i);
            labels.add(code.matches("[A-Za-z0-9_.\\-]+") ? code : "'" + code + "'");
        }
        final Set<String> seen = new HashSet<String>();
        final List<String> trees = new ArrayList<String>();
        final Map<String, Integer> qOf = new HashMap<String, Integer>();
        for (String binary : binaryTrees(labels)) {
            final Parameter keep = new Parameter.Default(1);
            final CharacterTreeModel model = new CharacterTreeModel("characterTree", binary, dataType, true, keep);
            final int internal = n - 1;
            for (int subset = 0; subset < (1 << internal); subset++) {
                for (int i = 0; i < internal; i++) {
                    keep.setParameterValueQuietly(i, (subset & (1 << i)) != 0 ? 0.0 : 1.0);
                }
                keep.setParameterValueQuietly(model.getSkeleton().getRoot().getNumber() - n, 1.0);
                final String contracted = model.getContraction().contractedNewick(null);
                if (seen.add(contracted)) {
                    trees.add(contracted);
                    qOf.put(contracted, model.getContraction().getContractedInternalNodeCount() - 1);
                }
            }
        }
        Collections.sort(trees, new Comparator<String>() {
            public int compare(String a, String b) {
                final int byQ = Integer.compare(qOf.get(a), qOf.get(b));
                return byQ != 0 ? byQ : a.compareTo(b);
            }
        });
        return trees;
    }

    private static List<String> binaryTrees(List<String> labels) {
        List<String> trees = new ArrayList<String>();
        trees.add("(" + labels.get(0) + "," + labels.get(1) + ")");
        for (int k = 2; k < labels.size(); k++) {
            final List<String> next = new ArrayList<String>();
            for (String t : trees) {
                next.addAll(insertEverywhere(t, labels.get(k)));
            }
            trees = next;
        }
        final List<String> out = new ArrayList<String>();
        for (String t : trees) {
            out.add(t + ";");
        }
        return out;
    }

    private static List<String> insertEverywhere(String tree, String tip) {
        final List<String> out = new ArrayList<String>();
        out.add("(" + tree + "," + tip + ")");
        if (tree.startsWith("(")) {
            int depth = 0;
            for (int i = 1; i < tree.length() - 1; i++) {
                final char c = tree.charAt(i);
                if (c == '(') depth++;
                if (c == ')') depth--;
                if (c == ',' && depth == 0) {
                    final String left = tree.substring(1, i), right = tree.substring(i + 1, tree.length() - 1);
                    for (String l : insertEverywhere(left, tip)) {
                        out.add("(" + l + "," + right + ")");
                    }
                    for (String r : insertEverywhere(right, tip)) {
                        out.add("(" + left + "," + r + ")");
                    }
                    break;
                }
            }
        }
        return out;
    }

    public List<String> getTreeList() {
        return treeList;
    }

    private long fillMasks(NodeRef node, List<Long> internalMasks) {
        if (tree.isExternal(node)) {
            return 1L << AbstractHierarchicalRateProvider.stateOfCode(tree.getDataType(),
                    tree.getNodeTaxon(node).getId());
        }
        long mask = 0;
        for (int c = 0; c < tree.getChildCount(node); c++) {
            mask |= fillMasks(tree.getChild(node, c), internalMasks);
        }
        if (!tree.isRoot(node)) {
            internalMasks.add(mask);
        }
        return mask;
    }

    private void compute() {
        if (known) {
            return;
        }
        Arrays.fill(values, 0.0);
        final List<Long> internalMasks = new ArrayList<Long>();
        fillMasks(tree.getRoot(), internalMasks);
        int k = 0;
        values[k++] = internalMasks.size();
        values[k++] = internalMasks.isEmpty() ? 1.0 : 0.0;
        if (indexOfTree != null) {
            final Integer index = indexOfTree.get(tree.getContractedNewick());
            values[k++] = index == null ? -1 : index;
        }
        for (long mask : internalMasks) {
            final Integer column = columnOfClade.get(mask);
            if (column != null) {
                values[column] = 1.0;
            }
        }
        k += cladeMasks.length;
        if (exchangeabilities) {
            model.getInfinitesimalMatrix(q);
            final double[] pi = model.getFrequencyModel().getFrequencies();
            for (int i = 0; i < stateCount; i++) {
                for (int j = i + 1; j < stateCount; j++) {
                    values[k++] = q[i * stateCount + j] / pi[j];
                }
            }
        }
        known = true;
    }

    @Override
    public int getDimension() {
        return values.length;
    }

    @Override
    public String getDimensionName(int dim) {
        return names.get(dim);
    }

    @Override
    public double getStatisticValue(int dim) {
        if (dim == 0) {
            known = false;
        }
        compute();
        return values[dim];
    }
}
