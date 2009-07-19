package dr.inference.trace;

/**
 * @author Alexei Drummond
 */
public abstract class AbstractTraceList implements TraceList {
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
        double[] values = new double[getStateCount()];
        int offset = (getBurnIn() / getStepSize());

        if (traceStatistics == null) {
            traceStatistics = new TraceCorrelation[getTraceCount()];
        }

        Trace trace = getTrace(index);
        trace.getValues(offset, values);
        traceStatistics[index] = new TraceCorrelation(values, getStepSize());
    }

    public void setBurnIn(int burnIn) {
        traceStatistics = null;
    }

    abstract Trace getTrace(int index);

    private TraceCorrelation[] traceStatistics = null;
}
