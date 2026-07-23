/*
 * Filter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.trace;

/**
 * @author Walter Xie
 */
public class Filter {

    final protected TraceType traceType; // for consistency matter
    protected String[] in; // bound of double or integer filter, values of string filter

    public Filter(String[] in, TraceType traceType) {
        if (traceType.isNumber() && in.length != 2)
            throw new IllegalArgumentException("Numeric filter must have both lower and upper ! trace type = " + traceType);
        setIn(in);
        this.traceType = traceType;
    }

    public boolean isIn(Object value) {
        if (traceType.isNumber()) {
            // double or integer
            return ( (Double)value >= Double.parseDouble(in[0]) && (Double)value <= Double.parseDouble(in[1]));
        }
        // String
        for (Object t : in) {
            if (t.toString().equals(value.toString())) {
                return true;
            }
        }
        return false;
    }

    /**
     * either lower and upper limits of double or integer filter, or values of string filter
     * @return String[]
     */
    public String[] getIn() {
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
