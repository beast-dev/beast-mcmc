package dr.inference.trace;

/**
 * @author Walter Xie
 */
public class Filter {

//    public boolean[] selected; // a mark, length = trace.values.size() = valueCount, all must be true initially
    final protected TraceFactory.TraceType traceType; // for consistency matter
    protected String[] in; // bound of double or integer filter, values of string filter

    public Filter(String[] in, TraceFactory.TraceType traceType) {
        if (traceType != TraceFactory.TraceType.STRING && in.length != 2)
            throw new IllegalArgumentException("Double or integer filter should have 2 bounds ! trace type = " + traceType);
        setIn(in);
        this.traceType = traceType;
    }

    public boolean isIn(Object value) {
        if (traceType == TraceFactory.TraceType.DOUBLE) {
            // double or integer
            return ( (Double)value >= Double.parseDouble(in[0]) && (Double)value <= Double.parseDouble(in[1]));
         }
//        else if (traceType == TraceFactory.TraceType.INTEGER) {
//            return ( (Integer)value >= Integer.parseInt(in[0]) && (Integer)value <= Integer.parseInt(in[1]));
//        }
        for (Object t : in) {
            if (t.toString().equals(value.toString())) {
                return true;
            }
        }
        return false;
    }

    public String[] getIn() {
//        String[] inString = new String[in.length];
//        for (int i = 0; i < in.length; i++) {
//            inString[i] = in[i].toString();
//        }
        return this.in;
    }

    public void setIn(String[] in) {
        this.in = in;
    }

    public String getStatusMessage() {
        String message = /*traceName +*/ " is filtered";
//        if (traceType == TraceFactory.TraceType.DOUBLE) {
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
