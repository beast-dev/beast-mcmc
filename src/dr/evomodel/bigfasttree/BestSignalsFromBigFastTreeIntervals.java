/*
 * BestSignalsFromBigFastTreeIntervals.java
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

package dr.evomodel.bigfasttree;

import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public class BestSignalsFromBigFastTreeIntervals extends BigFastTreeIntervals {

    private final Events originalEvents;

    public BestSignalsFromBigFastTreeIntervals(TreeModel tree) {
        super("signalsFromBigFastIntervals", tree);
        originalEvents = new Events(tree.getNodeCount());
    }

    public BestSignalsFromBigFastTreeIntervals(String name, TreeModel tree) {
        super(name, tree);
        originalEvents = new Events(tree.getNodeCount());
    }

    @Override
    public void calculateIntervals() {

        boolean changed = !intervalsKnown;

        if (!intervalsKnown && (originalEvents != null)) {
            originalEvents.copyEvents(events);
        }

        super.calculateIntervals();

        if (changed && (originalEvents != null)) {
            // Find first affected interval
            /*int firstEvent = findFirstEvent(originalEvents, events);
            if (firstEvent < events.size()) {
                fireModelChanged(new IntervalChangedEvent.FirstAffectedInterval(firstEvent - 1));
            }*/

            //using the code below to return range of intervals that are changed
            //find first and last interval affected
            int[] changedEvents = findChangedEvents(originalEvents, events);
            if (changedEvents[0] < events.size() && changedEvents[1] < events.size()) {
                fireModelChanged(new IntervalChangedEvent.AffectedIntervals(changedEvents));
            }
        }
    }

    private int findFirstEvent(Events lhs, Events rhs) {
        int i = 0;
        for (; i < lhs.size(); ++i) {
            if (lhs.getNode(i) != rhs.getNode(i) || lhs.getInterval(i) != rhs.getInterval(i)) break;
        }
        return i;
    }

    private int[] findChangedEvents(Events lhs, Events rhs) {
        int[] changes = new int[2];
        changes[0] = -1;
        changes[1] = -1;
        for (int i = 0; i < lhs.size(); ++i) {
            if (lhs.getNode(i) != rhs.getNode(i) || lhs.getInterval(i) != rhs.getInterval(i)) {
                if (changes[0] == -1) {
                    changes[0] = i-1;
                } else {
                    changes[1] = i-1;
                }
            }
        }
        return changes;
    }
}
