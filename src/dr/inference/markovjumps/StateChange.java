/*
 * StateChange.java
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

package dr.inference.markovjumps;

/**
 * A class to represent a state-change in a continuous-time Markov chain.
 * This class is equivalent to ColourChange in terms of parameters, but labels events
 * more generically as either forward or backwards (instead of just backwards) in time.
 * <p/>
 * TODO Remove dr.evolution.tree.ColourChange (wrong package)
 *
 * @author Alexei Drummond
 * @author Marc A. Suchard
 *         <p/>
 *         This work is supported by NSF grant 0856099
 */

public class StateChange {

    // time of the change in colour
    private double time;

    // the state associated with this time
    private int state;

    private int previousState;

    public StateChange(StateChange change) {
        this(change.time, change.state);
    }

    public StateChange(double time, int state, int previousState) {
        this.time = time;
        this.state = state;
        this.previousState = previousState;
    }

    public StateChange(double time, int state) {
        this.time = time;
        this.state = state;
        this.previousState = -1;
    }

    public StateChange clone() {
        return new StateChange(time, state, previousState);
    }

    /**
     * @return the time of this change in colour
     */
    public final double getTime() {
        return time;
    }

    /**
     * @return the state associated with the time
     */
    public final int getState() {
        return state;
    }

    public final int getPreviousState() {
        return previousState;
    }

    public final void setPreviousState(int state) {
        this.previousState = state;
    }

    /**
     * @param state sets the state of this event
     */
    public final void setState(int state) {
        this.state = state;
    }

    public final void setTime(double time) {
        this.time = time;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("{").append(getTime()).append(",").append(getState()).append("}");
        return sb.toString();
    }
}
