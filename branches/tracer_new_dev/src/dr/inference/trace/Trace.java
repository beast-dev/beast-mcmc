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
public class Trace<T> {

    public static final int INITIAL_SIZE = 1000;
    public static final int INCREMENT_SIZE = 1000;
    
    //    private TraceType traceType = TraceType.CONTINUOUS;
    protected T[] values = (T[]) new Object[INITIAL_SIZE];
    protected int valueCount = 0;
    protected final String name;

//    public Trace(String name) {
//        this.name = name;
//    }

    public Trace(String name, int initialSize, T initValue) {
        this.name = name;
        if (initialSize > 0) {
            this.values = (T[]) new Object[initialSize];
        }
        values[0] = initValue; // make getTraceType() working
    }

    public Trace(String name, T[] values) {
        this.name = name;
        this.values = (T[]) new Object[values.length];
        valueCount = values.length;
        System.arraycopy(values, 0, this.values, 0, values.length);
    }

    /**
     * @param value the valued to be added
     */
    public void add(T value) {
        if (valueCount == values.length) {
            T[] newValues = (T[]) new Object[valueCount + INCREMENT_SIZE];
            System.arraycopy(values, 0, newValues, 0, values.length);
            values = newValues;
        }

        values[valueCount] = value;
        valueCount++;
    }

    /**
     * @param values the values to be added
     */
    public void add(T[] values) {
        for (T value : values) {
            this.add(value);
        }
        valueCount += values.length;
    }

    public int getCount() {
        return valueCount;
    }

    public T getValue(int index) {
        return values[index];
    }

    public T[] createValues(int start, int length) {
        T[] destination = (T[]) new Object[length];
        getValues(start, destination, 0);
        return destination;
    }

    public void getValues(int start, T[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, valueCount - start);
    }

    public void getValues(int start, int count, T[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, count);
    }

    public String getName() {
        return name;
    }

//    public static Double[] arrayCopy(double[] src) {
//        Double[] dest = new Double[src.length];
//        for (int i=0; i < src.length; i++) {
//            dest[i] = Double.valueOf(src[i]);
//        }
//        return dest;
//    }

    public static double[] arrayConvert(Double[] src) {
        double[] dest = new double[src.length];
        for (int i=0; i < src.length; i++) {
            dest[i] = src[i].doubleValue();
        }
        return dest;
    }

    public static int[] arrayConvert(Integer[] src) {
        int[] dest = new int[src.length];
        for (int i=0; i < src.length; i++) {
            dest[i] = src[i].intValue();
        }
        return dest;
    }

    public static String[] arrayConvert(String[] src) {
        return src;
    }


    public Class getTraceType() {
        return values[0].getClass();
    }
//
//    public void setTraceType(TraceFactory.TraceType traceType) {
//        this.traceType = traceType;
//    }

}