package dr.inference.trace;

/**
 *
 */
public class Filter<T> {

    String traceName;
    T[] in; // it is selected, so that Trace.notSelected[i] = false

    public Filter(String traceName, T[] in) {
        this.traceName = traceName;
        this.in = in;
    }

    boolean isIn(T value) {

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

}
