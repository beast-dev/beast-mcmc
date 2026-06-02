/*
 * FilteredTraceList.java
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

package dr.inference.trace;

/**
 * not available for CombinedTraces yet
 */
public abstract class FilteredTraceList implements TraceList {

    // length = values[].size; the flag of which row is filtered,
    // where there is at least one column (parameter) having the value not in a selected range.
    // used in AbstractTraceList trace.getValues
    protected boolean[] filtered;

    protected void initFlag() {
        if (getTrace(0) != null) {
            filtered = new boolean[getTrace(0).getValueCount()];
        } else {
            throw new RuntimeException("Cannot initial filters ! getTrace(0) failed !");
        }
    }

    protected void updateFlag() {
        if (filtered == null) initFlag();

        for (int i=0; i < getTraceCount(); i++) {
            Trace trace = getTrace(i);
            Filter filter = trace.getFilter();
            if (filter != null) {
                if (trace.getValueCount() != filtered.length)
                    throw new RuntimeException("Invalid value size in Trace " + getTraceName(i));

                for (int j = 0; j < trace.getValueCount(); j++) {
                    // filter values not in
                    if ( ! filter.isIn(trace.getValue(j)) )
                        filtered[j] = true;
                }
            }
        }
    }

    public boolean hasAnyFilter() {
        for (int i=0; i < getTraceCount(); i++) {
            if (getTrace(i).getFilter() != null)
                return true;
        }
        return false;
    }

    public boolean hasFilter(int traceIndex) {
//        if (selected == null) return false;
        return getTrace(traceIndex).getFilter() != null;
    }

    public void setFilter(int traceIndex, Filter filter) {
        if (filtered == null) initFlag();
        getTrace(traceIndex).setFilter(filter);
        refreshStatistics();
    }

    public Filter getFilter(int traceIndex) {
        Filter filter = getTrace(traceIndex).getFilter();
        if (filtered == null && filter != null)
            throw new RuntimeException("The filter applied, but flag filtered[] is null ! " + getTraceName(traceIndex));
        return filter;
    }

    public void removeFilter(int traceIndex) {
        getTrace(traceIndex).setFilter(null);
        refreshStatistics();
    }

    public void removeAllFilters() {
        for (int i = 0; i < getTraceCount(); i++) {
            getTrace(i).setFilter(null);
        }
        filtered = null; // clean the flag
        refreshStatistics();
    }

    protected void refreshStatistics() {
        updateFlag();
        // must update filtered[] before analyseTrace
        for (int i = 0; i < getTraceCount(); i++) {
            analyseTrace(i);
        }
    }
}
