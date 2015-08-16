/*
 * FilteredTraceList.java
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

/**
 * not available for CombinedTraces yet
 */
public abstract class FilteredTraceList implements TraceList {

//    protected boolean[] selected; // length = values[].length = getValuesSize, all must be true initially

//    private void createSelected() { // will init in updateSelected()
//        if (getTrace(0) != null) {
//            selected = new boolean[getTrace(0).getValuesSize()];
//        } else {
//            throw new RuntimeException("Cannot initial filters ! getTrace(0) failed !");
//        }
//    }

//    private void initSelected() {
//        for (int i = 0; i < selected.length; i++) {
//            selected[i] = true;
//        }
//    }

    public boolean hasFilter(int traceIndex) {
//        if (selected == null) return false;
        return getTrace(traceIndex).getFilter() != null;
    }

    public void setFilter(int traceIndex, Filter filter) {
//        if (selected == null) createSelected();
        getTrace(traceIndex).setFilter(filter);
        refreshStatistics();
    }

    public Filter getFilter(int traceIndex) {
//        if (selected == null) return null;
        return getTrace(traceIndex).getFilter();
    }

    public void removeFilter(int traceIndex) {
        getTrace(traceIndex).setFilter(null);
        refreshStatistics();
    }

    public void removeAllFilters() {
        for (int i = 0; i < getTraceCount(); i++) {
            getTrace(i).setFilter(null);
        }
//        selected = null;
        refreshStatistics();// must be after "selected = null"
    }

    protected void refreshStatistics() {
//        updateSelected();
        for (int i = 0; i < getTraceCount(); i++) {
            analyseTrace(i);
        }
    }

//    private void updateSelected() {
//        if (selected != null) {
//            initSelected();
//            for (int traceIndex = 0; traceIndex < getTraceCount(); traceIndex++) {
//                if (getFilter(traceIndex) != null) {
//                    Trace trace = getTrace(traceIndex);
//                    if (trace.getValuesSize() != selected.length)
//                        throw new RuntimeException("updateSelected: length of values[] is different with selected[] in Trace "
//                                + getTraceName(traceIndex));
//
//                    for (int i = 0; i < trace.getValuesSize(); i++) {
//                        if (!trace.isIn(i)) { // not selected
//                            selected[i] = false;
//                        }
//                    }
//                }
//            }
//        }
//    }
}
