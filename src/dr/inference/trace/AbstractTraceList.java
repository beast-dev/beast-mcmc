package dr.inference.trace;

/**
 * @author Alexei Drummond
 */
public abstract class AbstractTraceList extends FilteredTraceList {
    public TraceDistribution getDistributionStatistics(int index) {
        return getCorrelationStatistics(index);
    }

    public TraceCorrelation getCorrelationStatistics(int index) {
        Trace trace = getTrace(index);
        if (trace == null) {
            return null;
        }
        return trace.getTraceStatistics();
    }

    public void analyseTrace(int index) {
        int start = (getBurnIn() / getStepSize());

//        if (traceStatistics == null) {
//            traceStatistics = new TraceCorrelation[getTraceCount()];
//            initFilters();
//        }

        Trace trace = getTrace(index);        
        TraceCorrelation traceCorrelation = new TraceCorrelation(
                trace.getValues(start, trace.getValuesSize(), selected),
                trace.getTraceType(), getStepSize());
        trace.setTraceStatistics(traceCorrelation);

//        System.out.println("index = " + index + " :  " + trace.getName() + "     " + trace.getTraceType());
    }

//    public void setBurnIn(int burnIn) {
//        traceStatistics = null;
//    }

//    abstract Trace getTrace(int index);

//    private TraceCorrelation[] traceStatistics = null;
}
