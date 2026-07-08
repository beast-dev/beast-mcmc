/*
 * MascotEventTape.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

/**
 * Builds a compact backward-time event list from a BEAST-X TreeModel and a
 * parameter of integer tip states ordered by external node index.
 */
public final class MascotEventTape {

    private final MascotCore.Event[] events;

    private MascotEventTape(MascotCore.Event[] events) {
        this.events = events;
    }

    public static MascotEventTape fromTree(TreeModel tree, Parameter tipStates, int stateCount) {
        int tipCount = tree.getExternalNodeCount();
        int internalCount = tree.getInternalNodeCount();
        if (tipStates.getDimension() != tipCount) {
            throw new IllegalArgumentException("tipStates dimension " + tipStates.getDimension() +
                    " does not match external node count " + tipCount);
        }

        MascotCore.Event[] events = new MascotCore.Event[tipCount + internalCount];
        int index = 0;
        for (int i = 0; i < tipCount; i++) {
            NodeRef node = tree.getExternalNode(i);
            int state = (int) Math.round(tipStates.getParameterValue(i));
            if (state < 0 || state >= stateCount) {
                throw new IllegalArgumentException("tip state out of range for external node " + i + ": " + state);
            }
            events[index++] = MascotCore.Event.sample(tree.getNodeHeight(node), node.getNumber(), state);
        }

        for (int i = 0; i < internalCount; i++) {
            NodeRef node = tree.getInternalNode(i);
            int childCount = tree.getChildCount(node);
            if (childCount != 2) {
                throw new IllegalArgumentException("MASCOT currently requires binary trees; node " +
                        node.getNumber() + " has " + childCount + " children");
            }
            NodeRef child1 = tree.getChild(node, 0);
            NodeRef child2 = tree.getChild(node, 1);
            events[index++] = MascotCore.Event.coalescence(
                    tree.getNodeHeight(node),
                    child1.getNumber(),
                    child2.getNumber(),
                    node.getNumber()
            );
        }

        return new MascotEventTape(events);
    }

    public MascotCore.Event[] getEvents() {
        return events.clone();
    }
}
