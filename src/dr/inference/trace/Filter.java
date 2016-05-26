/*
 * Filter.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
