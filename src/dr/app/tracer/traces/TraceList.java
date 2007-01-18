package dr.app.tracer.traces;

/**
 * An interface and default class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceList.java,v 1.2 2006/05/22 15:52:06 rambaut Exp $
 */

public interface TraceList {

	/** @return the name of this trace list */
	String getName();

	/** @return the number of traces in this trace list */
	int getTraceCount();

	/** @return the index of the trace with the given name
     *  @param name trace name
     */
	int getTraceIndex(String name);

	/** @return the name of the trace with the given index
     *  @param index trace ordinal index
     */
	String getTraceName(int index);

	/** @return the burn-in for this trace list */
	int getBurnIn();

	/** @return the number of states in the traces */
	int getStateCount();

	/** @return the size of the step between states */
	int getStepSize();

    /** @return the last state in the chain */
    int getMaxState();

    boolean isIncomplete();

    /** get the values of trace with the given index (without burnin) */
	void getValues(int index, double[] destination);

	/** get the values of trace with the given index (without burnin) */
	void getValues(int index, double[] destination, int offset);

	/** @return the trace distribution statistic object for the given index */
	TraceDistribution getDistributionStatistics(int index);

	/** @return the trace correlation statistic object for the given index */
	TraceCorrelation getCorrelationStatistics(int index);

	void analyseTrace(int index);
}