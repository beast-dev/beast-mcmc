package dr.inference.trace;

/**
 * not available for CombinedTraces yet
 */
public abstract class FilteredTraceList implements TraceList {

    protected boolean[] selected; // length = values[].length = valueCount, all must be true initially

    private void createSelected() { // will init in updateSelected()
        if (getTrace(0) != null) {
            selected = new boolean[getTrace(0).getValuesSize()];
        } else {
            throw new RuntimeException("Cannot initial filters ! getTrace(0) failed !");
        }
    }

    private void initSelected() {
        for (int i = 0; i < selected.length; i++) {
            selected[i] = true;
        }
    }

    public boolean hasFilter(int traceIndex) {
        if (selected == null) return false;
        return getTrace(traceIndex).getFilter() != null;
    }

    public void setFilter(int traceIndex, Filter filter) {
        if (selected == null) createSelected();
        getTrace(traceIndex).setFilter(filter);
        refreshStatistics();
    }

    public Filter getFilter(int traceIndex) {
        if (selected == null) return null;
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
        selected = null;
        refreshStatistics();// must be after "selected = null"
    }

    protected void refreshStatistics() {
        updateSelected();
        for (int i = 0; i < getTraceCount(); i++) {
            analyseTrace(i);
        }
    }

    private void updateSelected() {
        if (selected != null) {
            initSelected();
            for (int traceIndex = 0; traceIndex < getTraceCount(); traceIndex++) {
                if (getFilter(traceIndex) != null) {
                    Trace trace = getTrace(traceIndex);
                    if (trace.getValuesSize() != selected.length)
                        throw new RuntimeException("updateSelected: length of values[] is different with selected[] in Trace "
                                + getTraceName(traceIndex));

                    for (int i = 0; i < trace.getValuesSize(); i++) {
                        if (!trace.isIn(i)) { // not selected
                            selected[i] = false;
                        }
                    }
                }
            }
        }
    }
}
