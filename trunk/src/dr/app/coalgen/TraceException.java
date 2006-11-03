// TraceException.java

package dr.app.coalgen;

/**
 * An exception for traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceException.java,v 1.2 2005/06/15 17:20:54 rambaut Exp $
 */
 
public class TraceException extends Exception { 
	public TraceException() { super(); }
	public TraceException(String message) { super(message); }
}