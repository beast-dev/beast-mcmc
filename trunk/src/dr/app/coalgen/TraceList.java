// TraceList.java

package dr.app.coalgen;



/**
 * An interface and default class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceList.java,v 1.2 2004/10/01 22:40:02 alexei Exp $
 */
 
public interface TraceList {

	/** @return the name of this trace list */
	String getName();

	/** @return the number of traces in this trace list */
	int getTraceCount();
	
	/** @return the index of the trace with the given name */
	int getTraceIndex(String name);
	     
	/** @return the name of the trace with the given index */
	String getTraceName(int index);
	
	/** @return the burn-in for this trace list */
	int getBurnIn();

	/** @return the number of states in the traces */
	int getStateCount();
	
	/** @return the size of the step between states */
	int getStepSize();
	
	/** get the values of trace with the given index (without burnin) */
	void getValues(int index, double[] destination);
}