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

    public class WholeTree implements TreeChangedEvent {

        @Override public int getIndex() { return -1; }

        @Override public NodeRef getNode() { return null; }

        @Override public Parameter getParameter() { return null; }

        @Override public boolean isNodeChanged() { return false; }

        @Override public boolean isTreeChanged() { return true; }

        @Override public boolean isNodeParameterChanged() { return false; }

        @Override public boolean isHeightChanged() { return false; }
    }
}
