/*
 * AbstractTreeOperator.java
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

package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Statistic;
import dr.inference.operators.SimpleMCMCOperator;

import java.util.ArrayDeque;
import java.util.Deque;

/**
 * @author Andrew Rambaut
 */
public abstract class AbstractTreeOperator extends SimpleMCMCOperator implements Loggable {
    private final static int WINDOW_SIZE = 1000;
    private final Deque<Integer> topologyChanged = new ArrayDeque<>();
    private final Deque<Double> delta = new ArrayDeque<>();

    /* exchange sub-trees whose root are i and j */
    protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
                                 NodeRef iP, NodeRef jP) {

        tree.beginTreeEdit();
        tree.removeChild(iP, i);
        tree.removeChild(jP, j);
        tree.addChild(jP, i);
        tree.addChild(iP, j);

        tree.endTreeEdit();
    }

    /**
     * @param tree   the tree
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {
        if( tree.getChild(parent, 0) == child ) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
    }

    protected void setTopologyChanged(boolean changed) {
        setTopologyChanged(changed, 0.0);
    }

    protected void setTopologyChanged(boolean changed, double delta) {
        this.topologyChanged.addLast(changed ? 1 : 0);
        if (this.topologyChanged.size() > WINDOW_SIZE) {
            this.topologyChanged.removeFirst();
        }
        this.delta.addLast(delta);
        if (this.delta.size() > WINDOW_SIZE) {
            this.delta.removeFirst();
        }
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[] {
                new LogColumn.Abstract("topologyChanged") {
                    @Override
                    protected String getFormattedValue() {
                        int sum = 0;
                        for (int changed : topologyChanged) {
                            sum += changed;
                        }
                        return String.valueOf((double) sum / (double) (topologyChanged.size()));
                    }
                },
                new LogColumn.Abstract("mean_delta") {
                    @Override
                    protected String getFormattedValue() {
                        double sum = 0;
                        for (double delta : topologyChanged) {
                            sum += delta;
                        }
                        return String.valueOf(sum / (double) (topologyChanged.size()));
                    }
                }
        };

        return columns;
    }
}
