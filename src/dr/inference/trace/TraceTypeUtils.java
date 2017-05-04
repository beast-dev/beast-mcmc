package dr.inference.trace;

import java.util.List;

/**
 * @author Walter Xie
 */
public class TraceTypeUtils {

    /**
     * return true, if all selected traces are discrete
     * {@link dr.inference.trace.TraceType#isDiscrete() isDiscrete}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allDiscrete(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isDiscrete()) return false;
            }
        }
        return true;
    }

    /**
     * return true, if all selected traces are categorical
     * {@link dr.inference.trace.TraceType#isCategorical() isCategorical}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allCategorical(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isCategorical()) return false;
            }
        }
        return true;
    }

    public static boolean anyCategorical(List<TraceList> traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            if (traceNames == null) {
                for (int i = 0; i < traceList.getTraceCount(); i++) {
                    TraceDistribution td = traceList.getCorrelationStatistics(i);
                    if (td != null && ! td.getTraceType().isCategorical())
                        return true;
                }
            } else {
                for (String traceName : traceNames) {
                    int index = traceList.getTraceIndex(traceName);
                    TraceDistribution td = traceList.getCorrelationStatistics(index);
                    if (td != null && ! td.getTraceType().isCategorical())
                        return true;
                }
            }
        }
        return false;
    }


    /**
     * return true, if all selected traces are numbers
     * {@link dr.inference.trace.TraceType#isNumber() isNumber}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allNumeric(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isNumber()) return false;
            }
        }
        return true;
    }


    public static boolean anyNumeric(List<TraceList> traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            if (traceNames == null) {
                for (int i = 0; i < traceList.getTraceCount(); i++) {
                    TraceDistribution td = traceList.getCorrelationStatistics(i);
                    if (td != null && ! td.getTraceType().isNumber())
                        return true;
                }
            } else {
                for (String traceName : traceNames) {
                    int index = traceList.getTraceIndex(traceName);
                    TraceDistribution td = traceList.getCorrelationStatistics(index);
                    if (td != null && ! td.getTraceType().isNumber())
                        return true;
                }
            }
        }
        return false;
    }

}
