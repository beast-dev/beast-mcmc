package dr.app.tracer.traces;

/**
 * An exception for traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceException.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */

public class TraceException extends Exception {
	public TraceException() { super(); }
	public TraceException(String message) { super(message); }
}