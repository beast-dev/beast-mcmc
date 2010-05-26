package dr.inference.trace;

/**
 *
 */
public class Filter<T> {

    String traceName;
    TraceFactory.TraceType traceType;
    T[] in;

    boolean[] selected;

    public Filter(String traceName, TraceFactory.TraceType traceType, T[] in, T[] values) {
        this.traceName = traceName;
        this.traceType = traceType;
        this.in = in;

        createSelected(values);
    }

    public boolean isIn(T value) {

        if (value instanceof Double) {
            return ( (Double)value >= (Double)in[0] && (Double)value <= (Double)in[1]);
        }
        for (T t : in) {
            if (t.equals(value)) {
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
