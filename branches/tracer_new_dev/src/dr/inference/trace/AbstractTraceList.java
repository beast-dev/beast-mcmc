package dr.inference.trace;

/**
 * @author Alexei Drummond
 */
public abstract class AbstractTraceList extends FilteredTraceList {
    public TraceDistribution getDistributionStatistics(int index) {
        return getCorrelationStatistics(index);
    }

    public TraceCorrelation getCorrelationStatistics(int index) {
        if (traceStatistics == null) {
            return null;
        }

        return traceStatistics[index];
    }

    public void analyseTrace(int index) {
        int offset = (getBurnIn() / getStepSize());

        if (traceStatistics == null) {
            traceStatistics = new TraceCorrelation[getTraceCount()];
        }

        Trace trace = getTrace(index);        
        traceStatistics[index] = new TraceCorrelation(trace.createValues(offset, getStateCount()), getStepSize());
    }

    public void setBurnIn(int burnIn) {
        traceStatistics = null;
    }

//    abstract Trace getTrace(int index);

//    private TraceCorrelation[] traceStatistics = null;
}
