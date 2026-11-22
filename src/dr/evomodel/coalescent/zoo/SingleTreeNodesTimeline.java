/*
 * SingleTreeNodesTimeline.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.coalescent.zoo;
import dr.evomodel.coalescent.intervals.FastTreeIntervals;

// this builds a timeline with coalescent and sampling events for a single tree
public class SingleTreeNodesTimeline extends SingleTreeTimeline {
    private final int nNodes;

    private final double[] timeLine;
    private final boolean[] flagCoalescentEvent;
    private int[] numLineages;

    public SingleTreeNodesTimeline(FastTreeIntervals treeIntervals) {
        super(treeIntervals);
        this.nNodes = getnNodes();
        this.timeLine = new double[nNodes];
        this.numLineages = new int[nNodes];
        this.flagCoalescentEvent = new boolean[nNodes];
        makeLine();
    }

    protected void makeLine() {
        timeLine[0] = treeIntervals.getStartTime();
        flagCoalescentEvent[0] = false;
        for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
            timeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
            numLineages[nodeIndex] = treeIntervals.getLineageCount(nodeIndex);
            flagCoalescentEvent[nodeIndex] =
                    String.valueOf(treeIntervals.getIntervalType(nodeIndex - 1)).equals("coalescent"); //TODO this is hard coded ...
        }
    }


    protected int getnNodes() { return nNodes; }
    protected boolean[] getFlagCoalescentEvent() { return flagCoalescentEvent; }
    protected double[] getTimeLine() { return timeLine; }
}