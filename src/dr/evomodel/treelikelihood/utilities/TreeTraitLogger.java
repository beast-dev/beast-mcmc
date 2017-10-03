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

import java.util.ArrayList;
import java.util.List;

/**
 * A utility class to format traits on a tree in a LogColumn-based file.  This class takes an array of TraitProvider
 * and its package location is bound to move
 *
 * @author Marc A. Suchard
 */

public class TreeTraitLogger implements Loggable {

    public TreeTraitLogger(Tree tree,
                           TreeTrait[] traits) {
        this.tree = tree;
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

                for (int i = 0; i < tree.getNodeCount(); ++i) {
                    final NodeRef node = tree.getNode(i);

                    if (!tree.isRoot(node) || trait.getIntent() == TreeTrait.Intent.NODE) {

                        if (trait instanceof TreeTrait.DA) {

                            final int dim = ((double [])trait.getTrait(tree, node)).length;

                            for (int j = 0; j < dim; ++j) {

                                final int idx = j;
                                LogColumn column = new LogColumn.Abstract(trait.getTraitName() + "." + (i + 1) +
                                "." + (j + 1)) {
                                    @Override
                                    protected String getFormattedValue() {
                                        double[] x = (double[]) trait.getTrait(tree, node);
                                        return Double.toString(x[idx]);
                                    }
                                };

                                columns.add(column);
                            }

                        } else {

                            LogColumn column = new LogColumn.Abstract(trait.getTraitName() + "." + (i + 1)) {
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

    private Tree tree;
    private List<TreeTrait> treeTraits;
}
