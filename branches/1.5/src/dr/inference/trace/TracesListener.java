package dr.inference.trace;

/**
 * @author Alexei Drummond
 */
public interface TracesListener {

    void traceNames(String[] names);

    void traceRow(int state, double[] values);
}
