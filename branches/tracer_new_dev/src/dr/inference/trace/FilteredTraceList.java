package dr.inference.trace;

/**
 *
 */
public abstract class FilteredTraceList implements TraceList {

    Filter filter; // todo multi-filter for one traceList
    protected TraceCorrelation[] traceStatistics = null;

    public void setFilter(Filter filter) {
        this.filter = filter;
        traceStatistics[getTraceIndex(filter.getTraceName())].setFilter(filter);
    }

    public Filter getFilter() {
        return this.filter;
    }

    public Filter getFilter(String traceName) {
        if (traceStatistics != null && getTraceIndex(traceName) > 0) {
           this.filter = traceStatistics[getTraceIndex(traceName)].getFilter(); 
        } else {
           this.filter = null;
        }
//        doFilter();
        return this.filter;
    }

    public void removeFilter(String traceName) {
        this.filter = null;
        traceStatistics[getTraceIndex(traceName)].setFilter(null);
    }

    public void removeAllFilters() {
        this.filter = null;
        doFilter();
    }

    private void doFilter() {
        for (TraceDistribution traceD : traceStatistics) {
            traceD.setFilter(filter);
        }
    }

}
