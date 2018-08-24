/*
 * ArrayTraceList.java
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

package dr.inference.trace;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Alexei Drummond
 */
public class ArrayTraceList extends AbstractTraceList {

    String name;
    Map<String, Integer> traceIndex = new HashMap<String, Integer>();
    long burnin = 0;
    List<Trace> traces;
    long stepSize;

    public ArrayTraceList(String name, List<Trace> traces, int burnin) {
        this.name = name;
        this.traces = traces;
        this.burnin = burnin;

        for (int i = 0; i < traces.size(); i++) {
            traceIndex.put(traces.get(i).getName(), i);
        }

        Trace stateTrace = traces.get(0);
        this.stepSize = (int) Math.round((Double) stateTrace.getValue(1) - (Double) stateTrace.getValue(0));
    }

    public String getName() {
        return name;
    }

    public String getFullName() {
        return name;
    }
    
    public int getTraceCount() {
        return traces.size();
    }

    public int getTraceIndex(String name) {
        Integer index = traceIndex.get(name);
        if (index == null) return -1;
        return index;
    }

    public String getTraceName(int index) {
        return traces.get(index).getName();
    }

    /**
     * @return the burn-in for this trace list (the number of sampled states to discard)
     */
    public long getBurnIn() {
        return burnin;
    }

    /**
     * @return the number of states in the traces (after burnin removed)
     */
    public int getStateCount() {
        return traces.get(0).getValueCount();
    }

    /**
     * @return the number of states in the burnin
     */
    public int getBurninStateCount() {
        return (int) (getBurnIn() / (long) stepSize);
    }

    public long getStepSize() {
        return stepSize;
    }

    public long getMaxState() {
        return getStateCount() * getStepSize();
    }

    public boolean isIncomplete() {
        return false;
    }

    public List<Double> getValues(int index, int fromIndex, int toIndex) {
        throw new UnsupportedOperationException("not available");
    }

    public List<Double> getValues(int index) {
        throw new UnsupportedOperationException("not available");
    }

    public List<Double> getBurninValues(int index) {
        throw new UnsupportedOperationException("not available");
    }

    public Trace getTrace(int index) {
        return traces.get(index);
    }
}
