package dr.inference.trace;

/**
 *
 */
public class Filter<T> {

    Object[] in;

    public Filter(Object[] in) {
        setIn(in);
    }

    public boolean isIn(T value) {

        if (value instanceof Double) {
            // double filter passing String type for in[]
            return ( (Double)value >= Double.valueOf(in[0].toString()) && (Double)value <= Double.valueOf(in[1].toString()));
        }
        for (Object t : in) {
            if (t.toString().equals(value.toString())) {
                return true;
            }
        }
        return false;
    }
//
//    public String getTraceName() {
//        return traceName;
//    }
//
//    public TraceFactory.TraceType getTraceType() {
//        return traceType;
//    }

    public String[] getIn() {
        String[] inString = new String[in.length];
        for (int i = 0; i < in.length; i++) {
           inString[i] = in[i].toString();
        }
        return inString;
    }

    public void setIn(Object[] in) {
        this.in = in;
    }

    public String getStatusMessage() {
        String message = /*traceName +*/ " is filtered";
//        if (traceType == TraceFactory.TraceType.CONTINUOUS) {
//            message += " into range [" + in[0] + ", " + in[1] + "]";
//        } else {
//            message += " by selecting {";
//            for (Object t : in) {
//                message += t.toString() + ", ";
//            }
//            message = message.substring(0, message.lastIndexOf(", ")); // remove ", " for last in[]
//            message += "}";
//        }
        return message;
    }
}
