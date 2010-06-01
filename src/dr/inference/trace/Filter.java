package dr.inference.trace;

/**
 *
 */
public class Filter<T> {

    String traceName;
    TraceFactory.TraceType traceType;
    Object[] in;

    boolean[] selected;

    public Filter(String traceName, TraceFactory.TraceType traceType, Object[] in, T[] values) {
        this.traceName = traceName;
        this.traceType = traceType;

        setIn(in, values);
    }

    public boolean isIn(T value) {

        if (value instanceof Double) {
            return ( (Double)value >= (Double)in[0] && (Double)value <= (Double)in[1]);
        }
        for (Object t : in) {
            if (t.toString().equals(value.toString())) {
                return true;
            }
        }
        return false;
    }
    
    public String getTraceName() {
        return traceName;
    }

    public TraceFactory.TraceType getTraceType() {
        return traceType;
    }

    public String[] getIn() {
        String[] inString = new String[in.length];
        for (int i = 0; i < in.length; i++) {
           inString[i] = in[i].toString();
        }
        return inString;
    }

    public void setIn(Object[] in, T[] values) {
        this.in = in;
        createSelected(values);
    }

    public boolean getSelected(int index) {
        return selected[index];
    }


    private void createSelected(T[] values) {
        selected = new boolean[values.length]; // default = false

        for (int i = 0; i < values.length; i++) {
            if (isIn(values[i])) { // selected
                selected[i] = true;
            }
        }
    }

}
