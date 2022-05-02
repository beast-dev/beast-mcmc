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
