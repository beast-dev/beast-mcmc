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

import java.util.ArrayList;
import java.util.List;

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
        return values[index]; // filter?
    }

    public T[] getValues() {
        return values;
    }

    // used by others (e.g. CombinedTraces) without filter applied
    public void getValues(int start, T[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, valueCount - start);
    }

    // used by others (e.g. CombinedTraces) without filter applied
    public void getValues(int start, int count, T[] destination, int offset) {
        System.arraycopy(values, start, destination, offset, count);
    }

    public T[] getValues(int length, int start, int offset, boolean[] selected) {
        return this.getValues(length, start, offset, valueCount - start, selected);
    }

    public T[] getValues(int length, int start, int offset, int count, boolean[] selected) {
        T[] destination = (T[]) new Object[length];
        System.arraycopy(values, start, destination, offset, count);

        if (selected != null) {
            boolean[] destinationSelected = new boolean[length];
            System.arraycopy(selected, start, destinationSelected, offset, count);
            return getSeletedValues(destination, destinationSelected);
        } else {
            return destination;
        }
    }

    private T[] getSeletedValues(T[] values, boolean[] selected) {
        if (values.length != selected.length)
            throw new RuntimeException("getSeletedValues: length of values[] is different with selected[] in Trace " + name);

        List<T> valuesList = new ArrayList<T>();
        for (int i = 0; i < values.length; i++) {
            if (selected[i])
                valuesList.add(values[i]);
        }
        T[] afterSelected = (T[]) new Object[valuesList.size()];
        afterSelected = valuesList.toArray(afterSelected);
        return afterSelected;
    }

    public String getName() {
        return name;
    }

    public static <T> double[] arrayConvertToDouble(T[] src) {
        double[] dest = null;
        if (src != null) {
            dest = new double[src.length];
            for (int i = 0; i < dest.length; i++) {
                dest[i] = ((Number) src[i]).doubleValue();
            }
        }
        return dest;
    }

    public static <T> int[] arrayConvertToInt(T[] src) {
        int[] dest = null;
        if (src != null) {
            dest = new int[src.length];
            for (int i = 0; i < dest.length; i++) {
                dest[i] = ((Number) src[i]).intValue();
            }
        }
        return dest;
    }

    public static double[][] multiDArrayConvert(Double[][] src) {
        double[][] dest = null;
        if (src != null) {
            dest = new double[src.length][src[0].length];
            for (int i = 0; i < dest.length; i++) {
                for (int j = 0; j < dest[i].length; j++) {
                    dest[i][j] = src[i][j].doubleValue();
                }
            }
        }
        return dest;
    }

    public Class getTraceType() {
        return values[0].getClass();
    }
//
//    public void setTraceType(TraceFactory.TraceType traceType) {
//        this.traceType = traceType;
//    }

    //******************** TraceCorrelation ****************************
    protected TraceCorrelation<T> traceStatistics;

    public TraceCorrelation getTraceStatistics() {
        return traceStatistics;
    }

    public void setTraceStatistics(TraceCorrelation traceStatistics) {
        this.traceStatistics = traceStatistics;
    }

    //******************** Filter ****************************
    protected Filter<T> filter;

    public void setFilter(Filter filter) {
        this.filter = filter;
//        if (traceStatistics == null)
//            throw new RuntimeException("Cannot set filter because traceStatistics = null in Trace " + name);
//        traceStatistics =
    }

    public Filter getFilter() {
        return filter;
    }

    public boolean isIn(int i) {        
        return filter.isIn(values[i]);
    }
}