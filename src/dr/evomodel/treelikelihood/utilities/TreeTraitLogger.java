/*
 * TreeTraitLogger.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.xml.Report;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * A utility class to format traits on a tree in a LogColumn-based file.  This class takes an array of TraitProvider
 * and its package location is bound to move
 *
 * @author Marc A. Suchard
 */

public class TreeTraitLogger implements Loggable, Reportable {

    public TreeTraitLogger(Tree tree, TreeTrait[] traits) {
        this(tree, traits, NodeRestriction.ALL);
    }

    public TreeTraitLogger(Tree tree,
                           TreeTrait[] traits,
                           NodeRestriction nodeRestriction) {
        this.tree = tree;
        this.nodeRestriction = nodeRestriction;
        addTraits(traits);
    }

    private void addTraits(TreeTrait[] traits) {
        if (treeTraits == null) {
            treeTraits = new ArrayList<TreeTrait>();
        }

        for (TreeTrait trait : traits) {
            treeTraits.add(trait);
        }
    }

    public LogColumn[] getColumns() {

        List<LogColumn> columns = new ArrayList<LogColumn>();

        for (final TreeTrait trait : treeTraits) {

            if (trait.getIntent() == TreeTrait.Intent.WHOLE_TREE) {
                LogColumn column = new LogColumn.Abstract(trait.getTraitName()) {
                    @Override
                    protected String getFormattedValue() {
                        return trait.getTraitString(tree, null);
                    }
                };

                columns.add(column);

            } else {

                for (int i = nodeRestriction.begin(tree); i < nodeRestriction.end(tree); ++i) {
                    final NodeRef node = tree.getNode(i);

                    if (!tree.isRoot(node) || trait.getIntent() == TreeTrait.Intent.NODE) {

                        if (trait instanceof TreeTrait.DA) {

                            final int dim = ((double [])trait.getTrait(tree, node)).length;
                            for (int j = 0; j < dim; ++j) {

                                final int idx = j;
                                LogColumn column = new LogColumn.Abstract(trait.getTraitName()
                                        + "." + (i + 1) + "." + (j + 1)) {
                                    @Override
                                    protected String getFormattedValue() {
                                        double[] x = (double[]) trait.getTrait(tree, node);
                                        return Double.toString(x[idx]);
                                    }
                                };

                                columns.add(column);
                            }

                        } else {

                            LogColumn column = new LogColumn.Abstract(trait.getTraitName()
                                    + "." + (i + 1)) {
                                @Override
                                protected String getFormattedValue() {
                                    return trait.getTraitString(tree, node);
                                }
                            };

                            columns.add(column);
                        }
                    }
                }
            }
        }

        return columns.toArray(new LogColumn[columns.size()]);
    }

    private final Tree tree;
    private List<TreeTrait> treeTraits;
    private final NodeRestriction nodeRestriction;

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        for (LogColumn column : getColumns()) {
            sb.append(column.getFormatted());
            sb.append("\t");
        }

        return sb.toString();
    }

    public enum NodeRestriction {
        ALL {
            int begin(Tree tree) { return 0; }
            int end(Tree tree) { return tree.getNodeCount(); }
        },
        EXTERNAL {
            int begin(Tree tree) { return 0; }
            int end(Tree tree) { return tree.getExternalNodeCount(); }
        },
        INTERNAL {
            int begin(Tree tree) { return tree.getExternalNodeCount(); }
            int end(Tree tree) { return tree.getNodeCount(); }
        };

        abstract int begin(Tree tree);
        abstract int end(Tree tree);

        public static NodeRestriction parse(String text) {
             String lower = text.toLowerCase();
             if (lower.compareTo("external") == 0) {
                 return EXTERNAL;
             } else if (lower.compareTo("internal") == 0) {
                return INTERNAL;
            } else {
                 return ALL;
            }
        }
    }
}
