package dr.inference.trace;

/**
 *
 */
public abstract class FilteredTraceList implements TraceList {

    Filter[] filters; // store configured filter
    protected TraceCorrelation[] traceStatistics = null;
    int currentFilter = -1;

    public void setFilter(Filter filter) {
        int fId = getTraceIndex(filter.getTraceName());
        filters[fId] = filter;
        currentFilter = fId;
        doFilter(filter);

        createTraceFilter(filter);
    }

    public Filter getFilter() {
        if (currentFilter < 0) return null;
        return filters[currentFilter];
    }

    public Filter getFilter(String traceName) {
        int fId = getTraceIndex(traceName);
        if (traceStatistics != null && fId > 0) {
            return filters[fId];
        }
        return null;
    }

    public void removeFilter(String traceName) {
        int fId = getTraceIndex(traceName);
        if (filters[fId] != null) {
            filters[fId] = null;
            doFilter(null);
        }
        selected = null;
    }

    public void removeAllFilters() {
        for (Filter f : filters) {
            f = null;
        }
        doFilter(null);
        selected = null;
        currentFilter = -1;
    }

    private void doFilter(Filter filter) {
        for (TraceCorrelation traceC : traceStatistics) {
            traceC.setFilter(filter);
        }
    }

//    protected abstract void updateTraceList();

    protected void initFilters() { // used in void analyseTrace(int index)
        if (filters == null) filters = new Filter[traceStatistics.length]; // traceStatistics.length = getTraceCount()
    }


    //******************** Trace ****************************
    protected boolean[] selected; // length = valueCount

    public abstract void createTraceFilter(Filter filter);

    public void setTraceFilter(boolean[] selected) {
        this.selected = selected;
    }

    public boolean[] getTraceFilter() {
        return selected;
    }

}
