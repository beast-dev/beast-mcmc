package dr.inference.trace;

/**
 *
 */
public class FilteredTraceList implements TraceList {

    TraceList base;
    Filter filter; // todo Filter[]?
    TraceDistribution[] filteredTraces;

    public FilteredTraceList(TraceList baseList) {
        this.base = baseList;


    }

    public void setFilter(Filter filter) {
        this.filter = filter;
        doFilter();
    }

     public void removeFilter() {
        this.filter = null;
        doFilter();
    }

    private void doFilter() {

        for (TraceDistribution traceD : filteredTraces) {
            traceD.setFilter(filter);
        }
    }

    public String getName() {
        return base.getName();
    }

    public int getTraceCount() {
        return base.getTraceCount();
    }

    public int getTraceIndex(String name) {
        return base.getTraceIndex(name);
    }

    public String getTraceName(int index) {
        return base.getTraceName(index);
    }

    public int getBurnIn() {
        return base.getBurnIn();
    }

    public int getStateCount() {
        return base.getStateCount();
    }

    public int getBurninStateCount() {
        return base.getBurninStateCount();
    }

    public int getStepSize() {
        return base.getStepSize();
    }

    public int getMaxState() {
        return base.getMaxState();
    }

    public boolean isIncomplete() {
        return base.isIncomplete();
    }

    public <T> void getValues(int index, T[] destination) {
        base.getValues(index, destination);
    }

    public <T> void getValues(int index, T[] destination, int offset) {
        base.getValues(index, destination, offset);
    }

    public <T> void getBurninValues(int index, T[] destination) {
        base.getBurninValues(index, destination);
    }

    public TraceDistribution getDistributionStatistics(int traceIndex) {
        return base.getDistributionStatistics(traceIndex);
    }

    public TraceCorrelation getCorrelationStatistics(int traceIndex) {
        return base.getCorrelationStatistics(traceIndex);
    }

    public void analyseTrace(int index) {
        base.analyseTrace(index);
    }

    public Trace getTrace(int index) {
        return base.getTrace(index);
    }
}
