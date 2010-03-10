/*
 * Trace.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import java.lang.reflect.Array;


/**
 * A simple class that stores a trace for a single statistic
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Trace.java,v 1.11 2005/07/11 14:07:26 rambaut Exp $
 */
public abstract class Trace<T> implements TraceType<T> {

    //    private TraceType traceType = TraceType.CONTINUOUS;
    protected T[] values = (T[]) new Object[INITIAL_SIZE];
    protected int valueCount = 0;
    protected final String name;

//    public Trace(String name) {
//        this.name = name;
//    }

    public Trace(String name, int initialSize) {
        this.name = name;
        if (initialSize > 0) {
            this.values = (T[]) new Object[initialSize];
        }
    }

    public Trace(String name, T[] values) {
        this.name = name;
        this.values = (T[]) new Object[values.length];
        valueCount = values.length;
        System.arraycopy(values, 0, this.values, 0, values.length);
    }

    /**
     * @param values the values to be added
     */
//    public <T> void add(T[] values) {
//        for (T value : values) {
//            this.add(value);
//        }
//        valueCount += values.length;
//    }

    public int getCount() {
        return valueCount;
    }

    public T getValue(int index) {
        return values[index];
    }

    public void getValues(int start, double[] destination) {
        getValues(start, destination, 0);
    }

    public void getValues(int start, double[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, valueCount - start);
    }

    public void getValues(int start, int count, double[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, count);
    }

    public String getName() {
        return name;
    }    

//    public TraceType getTraceType() {
//        return traceType;
//    }
//
//    public void setTraceType(TraceType traceType) {
//        this.traceType = traceType;
//    }

//    public class DiscreteTrace extends Trace<Integer> {
//        public DiscreteTrace(String name, int initialSize) {
//            super(name, initialSize);
//        }
//
//        public DiscreteTrace(String name, Integer[] values) {
//            super(name, values);
//        }
//
//        public void add(Integer value) {
//            if (valueCount == values.length) {
//                Integer[] newValues = new Integer[valueCount + INCREMENT_SIZE];
//                System.arraycopy(values, 0, newValues, 0, values.length);
//                super.values = newValues;
//            }
//
//            super.values[valueCount] = value;
//            super.valueCount++;
//        }
//
//        public Integer parserValueWithType(String value) {
//            return Integer.parseInt(value);
//        }
//    }

//    public class ContinuousTrace extends Trace<Double> {
//        public ContinuousTrace(String name, int initialSize) {
//            super(name, initialSize);
//        }
//
//        public ContinuousTrace(String name, Double[] values) {
//            super(name, values);
//        }
//
//        public void add(Double value) {
//            if (valueCount == values.length) {
//                Double[] newValues = new Double[valueCount + INCREMENT_SIZE];
//                System.arraycopy(values, 0, newValues, 0, values.length);
//                super.values = newValues;
//            }
//
//            super.values[valueCount] = value;
//            super.valueCount++;
//        }
//
//        public Double parserValueWithType(String value) {
//            return Double.parseDouble(value);
//        }
//    }

//    public class CategoryTrace extends Trace<String> {
//        public CategoryTrace(String name, int initialSize) {
//            super(name, initialSize);
//        }
//
//        public CategoryTrace(String name, String[] values) {
//            super(name, values);
//        }
//
//        public void add(String value) {
//            if (valueCount == values.length) {
//                String[] newValues = new String[valueCount + INCREMENT_SIZE];
//                System.arraycopy(values, 0, newValues, 0, values.length);
//                super.values = newValues;
//            }
//
//            super.values[valueCount] = value;
//            super.valueCount++;
//        }
//
//        public String parserValueWithType(String value) {
//            return value;
//        }
//    }

}