/*
 * TreeChangedEvent.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.inference.model.Parameter;

/**
 * @author Alexei J. Drummond
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public interface TreeChangedEvent {
    int getIndex();

    NodeRef getNode();

    Parameter getParameter();

    boolean isNodeChanged();

    boolean isTreeChanged();

    boolean isNodeOrderChanged();

    boolean isNodeParameterChanged();

    boolean isHeightChanged();

    boolean isOnlyHeightChanged();

    class WholeTree implements TreeChangedEvent {

        @Override public int getIndex() { return -1; }

        @Override public NodeRef getNode() { return null; }

        @Override public Parameter getParameter() { return null; }

        @Override public boolean isNodeChanged() { return false; }

        @Override public boolean isNodeOrderChanged() { return true; } // TODO MAS: guessing on this

        @Override public boolean isTreeChanged() { return true; }

        @Override public boolean isNodeParameterChanged() { return false; }

        @Override public boolean isHeightChanged() { return false; }

        @Override public boolean isOnlyHeightChanged() { return false; }
    }

    class NodeOnTree implements TreeChangedEvent {

        private final NodeRef node;

        NodeOnTree(NodeRef node) {
            this.node = node;
        }

        @Override
        public int getIndex() { return -1; }

        @Override
        public NodeRef getNode() { return node; }

        @Override
        public Parameter getParameter() { return null; }

        @Override
        public boolean isNodeChanged() { return true; }

        @Override
        public boolean isNodeOrderChanged() { return true; } // TODO MAS: guessing on this

        @Override
        public boolean isTreeChanged() { return false; }

        @Override
        public boolean isNodeParameterChanged() { return false; }

        @Override
        public boolean isHeightChanged() { return false; }

        @Override public boolean isOnlyHeightChanged() { return false; }
    }

    static TreeChangedEvent create() {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return null; }

            @Override public Parameter getParameter() { return null; }

            @Override public boolean isNodeOrderChanged() { return true; }

            @Override public boolean isNodeChanged() { return false; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return false; }

            @Override public boolean isHeightChanged() { return false; }

            @Override public boolean isOnlyHeightChanged() { return false; }
        };
    }

    static TreeChangedEvent create(final boolean isNodeOrderChanged, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return null; }

            @Override public Parameter getParameter() { return null; }

            @Override public boolean isNodeChanged() { return false; }

            @Override public boolean isNodeOrderChanged() { return isNodeOrderChanged; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return false; }

            @Override public boolean isHeightChanged() { return isHeightChanged; }

            @Override public boolean isOnlyHeightChanged() { return false; }
        };
    }


    static TreeChangedEvent create(final NodeRef node, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return null; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeOrderChanged() { return isHeightChanged; }

            @Override public boolean isNodeParameterChanged() { return false; }

            @Override public boolean isHeightChanged() { return isHeightChanged; }

            @Override public boolean isOnlyHeightChanged() { return false; }
        };
    }

    static TreeChangedEvent create(final NodeRef node, final Parameter parameter, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return parameter; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeOrderChanged() { return isHeightChanged; }

            @Override public boolean isNodeParameterChanged() { return true; }

            @Override public boolean isHeightChanged() { return isHeightChanged; }

            @Override public boolean isOnlyHeightChanged() { return false; }
        };
    }

    static TreeChangedEvent create(final NodeRef node, final Parameter parameter, final int index, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return index; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return parameter; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeOrderChanged() { return isHeightChanged; }

            @Override public boolean isNodeParameterChanged() { return true; }

            @Override public boolean isHeightChanged() { return isHeightChanged; }

            @Override public boolean isOnlyHeightChanged() { return isHeightChanged; }
        };
    }
}
