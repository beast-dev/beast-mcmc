/*
 * TraceFactory.java
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
 * @author Alexei Drummond
 * @author Walter Xie
 * @author Andrew Rambaut
 */
public class TraceFactory {

    public enum TraceType {
        // changed this to 'real' as this is less Comp Sci. than 'double'
        DOUBLE("real", "R", Double.class),
        INTEGER("integer", "I", Integer.class),
        STRING("categorical", "C", String.class);

        TraceType(String name, String brief, Class type) {
            this.name = name;
            this.brief = brief;
            this.type = type;
        }

        public String toString() {
            return name;
        }

        public String getBrief() {
            return brief;
        }

        public Class getType() {
            return type;
        }

        private final String name;
        private final String brief;
        private final Class type;
    }

//    public static Trace createTrace(TraceType traceType, String name, int initialSize) {
//
////        Double[] d = new Double[10];
////        Double[] t = new Double[10];
////        System.arraycopy(d, 0, t, 0, d.length);
//
//        // System.out.println("create trace (" + name + ") with type " + traceType);
//
//        switch (traceType) {
//            case DOUBLE:
//                return new Trace<Double>(name, initialSize, (double) 0);
//            case INTEGER:
//                return new Trace<Integer>(name, initialSize, 0);
//            case STRING:
//                return new Trace<String>(name, initialSize, "initial_value");
//        }
//        throw new IllegalArgumentException("The trace type " + traceType + " is not recognized.");
//    }

}

