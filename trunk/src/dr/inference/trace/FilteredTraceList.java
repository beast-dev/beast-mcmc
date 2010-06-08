package dr.inference.trace;

/**
 *
 */
public abstract class FilteredTraceList implements TraceList {

    Filter[] filters; // store configured filter
    protected TraceCorrelation[] traceStatistics = null;

    public void setFilter(Filter filter) {
        int fId = getTraceIndex(filter.getTraceName());
        filters[fId] = filter;
        doFilter(filter);
    }

    public Filter getFilter(int index) {
        return filters[index];
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
    }

    public void removeAllFilters() {
        for (Filter f : filters) {
            f = null;
        }
        doFilter(null);
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
