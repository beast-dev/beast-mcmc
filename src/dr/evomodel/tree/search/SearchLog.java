/*
 * SearchLog.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * The search's outputs (plan 7 §12): a tab-separated log with one row per phase or round, and the final tree as
 * NEXUS, with each node's height (and branch rate) and F with its terms in the tree's comment.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SearchLog {

    public static final String[] COLUMNS = {"start", "round", "phase", "nT0", "nT1", "nT2", "nAccepted", "F_before",
            "F_after", "F_best", "move", "dF0", "dF1", "dF2", "rootedRF_to_start", "topology", "evaluations",
            "maximizer_reason", "guard_reverts", "seconds"};

    private final TreeModel tree;
    private final ArbitraryBranchRates rates;
    private final Likelihood objective;
    private final List<String> rows = new ArrayList<String>();

    SearchLog(TreeModel tree, ArbitraryBranchRates rates, Likelihood objective) {
        this.tree = tree;
        this.rates = rates;
        this.objective = objective;
    }

    void clear() {
        rows.clear();
    }

    void row(int round, String phase, int nT0, int nT1, int nT2, int nAccepted, double fBefore, double fAfter,
             double fBest, String move, double dF0, double dF1, double dF2, Tree startTree, double seconds,
             long evaluations, String maximizerReason, int guardReverts) {
        rows.add(join(1, round, phase, nT0, nT1, nT2, nAccepted, number(fBefore), number(fAfter), number(fBest),
                move == null ? "" : move, number(dF0), number(dF1), number(dF2),
                SplitDistance.rooted(startTree, tree), TopologyHash.toHex(TopologyHash.of(tree)), evaluations,
                maximizerReason, guardReverts, String.format(Locale.US, "%.3f", seconds)));
    }

    /**
     * @return the rows so far, tab separated, without the header
     */
    public List<String> getRows() {
        return new ArrayList<String>(rows);
    }

    void writeLog(File file) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println(join((Object[]) COLUMNS));
        for (String row : rows) {
            writer.println(row);
        }
        writer.close();
    }

    void writeTree(File file, double f, String stopReason) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println("#NEXUS");
        writer.println();
        writer.println("Begin taxa;");
        writer.println("\tDimensions ntax=" + tree.getTaxonCount() + ";");
        writer.println("\tTaxlabels");
        for (int i = 0; i < tree.getTaxonCount(); ++i) {
            writer.println("\t\t" + quote(tree.getTaxonId(i)));
        }
        writer.println("\t\t;");
        writer.println("End;");
        writer.println();
        writer.println("Begin trees;");
        writer.println("\tTranslate");
        for (int i = 0; i < tree.getTaxonCount(); ++i) {
            writer.println("\t\t" + (i + 1) + " " + quote(tree.getTaxonId(i)) +
                    (i < tree.getTaxonCount() - 1 ? "," : ""));
        }
        writer.println("\t\t;");
        StringBuilder comment = new StringBuilder("&F=").append(number(f));
        if (objective instanceof CompoundLikelihood) {
            for (Likelihood term : ((CompoundLikelihood) objective).getLikelihoods()) {
                if (term.getId() != null) {
                    comment.append(",").append(term.getId()).append("=").append(number(term.getLogLikelihood()));
                }
            }
        }
        comment.append(",objective=").append(objective.getId() != null ? objective.getId() : "objective");
        comment.append(",stop=\"").append(stopReason).append("\"");
        StringBuilder newick = new StringBuilder();
        newick(tree.getRoot(), newick);
        writer.println("tree TREE1 [" + comment + "] = [&R] " + newick + ";");
        writer.println("End;");
        writer.close();
    }

    private void newick(NodeRef node, StringBuilder sb) {
        if (tree.isExternal(node)) {
            sb.append(tree.getTaxonIndex(tree.getNodeTaxon(node).getId()) + 1);
        } else {
            sb.append("(");
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                if (c > 0) {
                    sb.append(",");
                }
                newick(tree.getChild(node, c), sb);
            }
            sb.append(")");
        }
        sb.append("[&height=").append(number(tree.getNodeHeight(node)));
        if (rates != null && !tree.isRoot(node)) {
            sb.append(",rate=").append(number(rates.getBranchRate(tree, node)));
        }
        sb.append("]");
        if (!tree.isRoot(node)) {
            sb.append(":").append(number(tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node)));
        }
    }

    private static String quote(String id) {
        return id.matches("[A-Za-z0-9_.]+") ? id : "'" + id.replace("'", "''") + "'";
    }

    private static String number(double value) {
        return Double.isNaN(value) ? "" : String.valueOf(value);
    }

    private static String join(Object... values) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.length; ++i) {
            if (i > 0) {
                sb.append('\t');
            }
            sb.append(values[i]);
        }
        return sb.toString();
    }
}
