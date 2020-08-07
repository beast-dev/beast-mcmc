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

    boolean isNodeParameterChanged();

    boolean isHeightChanged();

    static TreeChangedEvent create() {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return null; }

            @Override public Parameter getParameter() { return null; }

            @Override public boolean isNodeChanged() { return false; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return false; }

            @Override public boolean isHeightChanged() { return false; }
        };
    }
    public static TreeChangedEvent create(final NodeRef node, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return null; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return false; }

            @Override public boolean isHeightChanged() { return isHeightChanged; }
        };
    }

    public static TreeChangedEvent create(final NodeRef node, final Parameter parameter, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return -1; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return parameter; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return true; }

            @Override public boolean isHeightChanged() { return false; }
        };
    }

    public static TreeChangedEvent create(final NodeRef node, final Parameter parameter, final int index, final boolean isHeightChanged) {
        return new TreeChangedEvent() {
            @Override public int getIndex() { return index; }

            @Override public NodeRef getNode() { return node; }

            @Override public Parameter getParameter() { return parameter; }

            @Override public boolean isNodeChanged() { return true; }

            @Override public boolean isTreeChanged() { return true; }

            @Override public boolean isNodeParameterChanged() { return true; }

            @Override public boolean isHeightChanged() { return false; }
        };
    }

}
