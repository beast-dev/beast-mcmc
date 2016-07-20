/*
 * Trace.java
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

import java.util.*;

/**
 * A simple class that stores a trace for a single statistic
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Trace.java,v 1.11 2005/07/11 14:07:26 rambaut Exp $
 */
public class Trace<T> { // TODO get rid of generic to make things easy

//    public static final int INITIAL_SIZE = 1000;
//    public static final int INCREMENT_SIZE = 1000;

    // use <Double> for integer, but traceType must = ORDINAL, because of legacy issue at analyseCorrelationContinuous
    protected TraceType traceType = TraceType.REAL;
    protected List<T> values = new ArrayList<T>(); // TODO change to String only, and parse to double, int or string in getValues according to trace type
    //    protected int valueCount = 0;
    protected String name;

//    private Object[] range;

    public Trace(String name) { // traceType = TraceFactory.TraceType.DOUBLE; 
        this.name = name;
    }

    public Trace(String name, TraceType traceType) {
        this.name = name;
        setTraceType(traceType);
    }

//    public Trace(String name, T[] valuesArray) {
//        this(name);
////        List<T> newVL = Arrays.asList(valuesArray);
//        Collections.addAll(this.values, valuesArray);
//    }

    /**
     * @param value the valued to be added
     */
    public void add(T value) {
        values.add(value);
    }

    /**
     * @param valuesArray the values to be added
     */
    public void add(T[] valuesArray) {
        Collections.addAll(this.values, valuesArray);
    }

    public int getValueCount() {
        return values.size();
    }

    public T getValue(int index) {
//        TODO filtered ?
//        if (getFilter() != null && !getFilter().isIn(index)) {
//           return null; // filtered
//        }
        return values.get(index);
    }

    public double[] getRange() { // Double => bounds; Integer and String => unique values

        if (getValueCount() < 1) throw new IllegalArgumentException("Cannot find values in trace " + getName());

        if (getTraceType().isNumber()) {

            Double min = Double.MAX_VALUE;
            Double max = Double.MIN_VALUE;
            for (Object t : values) {
                if ((Double) t < min) {
                    min = (Double) t;
                } else if ((Double) t > max) {
                    max = (Double) t;
                }
            }
            return new double[] {min, max};
        } else {
            throw new UnsupportedOperationException("cannot call getRange for categorical data");
        }
    }

    public Set<String> getCategoricalValues() {
        return new TreeSet<String>((List<String>) values);
    }
    /**
     * @param fromIndex low endpoint (inclusive) of the subList.
     * @param toIndex   high endpoint (exclusive) of the subList.
     * @return The list of values (which are selected values if filter applied)
     */
    public List<T> getValues(int fromIndex, int toIndex) {
        if (toIndex > getValueCount() || fromIndex > toIndex)
            throw new RuntimeException("Invalid index : fromIndex = " + fromIndex + "; toIndex = " + toIndex
                    + "; List size = " + getValueCount() + "; in Trace " + name);

        if (getFilter() == null) {
            return values.subList(fromIndex, toIndex);
        } else {
//            if (filter.selected.length != getValuesSize())
//                throw new IllegalArgumentException("Trace " + name + " size of values is different with filter selected[] ! ");

            List<T> valuesList = new ArrayList<T>();
            for (int i = fromIndex; i < toIndex; i++) {
                if (getFilter().isIn(values.get(i)))
                    valuesList.add(values.get(i));
            }

            if (valuesList.size() < 1) throw new RuntimeException("There is no value left after applying filter !");

            return valuesList;
        }
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

//    public Class getTraceType() {
//        if (values.get(0) == null) {
//            return null;
//        }
//        return values.get(0).getClass();
//    }

    public TraceType getTraceType() {
        return traceType;
    }

    public void setTraceType(TraceType traceType) {
        this.traceType = traceType;
    }

    //******************** TraceCorrelation ****************************
    protected TraceCorrelation<T> traceStatistics;

    public TraceCorrelation getTraceStatistics() {
        return traceStatistics;
    }

    public void setTraceStatistics(TraceCorrelation traceStatistics) {
        this.traceStatistics = traceStatistics;
    }

    //******************** Filter ****************************
    protected Filter filter = null;

    public void setFilter(Filter filter) {
        this.filter = filter;
    }

    public Filter getFilter() {
        return filter;
    }

    //******************** Trace Double ****************************
/*    public class D extends Trace<Double> {

        public D(String name, int initialSize) {
            super.name = name;
            super.values = new Double[initialSize];
//            values[0] = initValue; // make getTraceType() working
        }

        public D(String name, Double[] values) {
            this(name, values.length);
            valueCount = values.length;
            System.arraycopy(values, 0, this.values, 0, values.length);
        }

        public Double[] getValues(int length, int start, int offset, boolean[] selected) {
            return this.getValues(length, start, offset, valueCount - start, selected);
        }

        public Double[] getValues(int length, int start, int offset, int count, boolean[] selected) {
            Double[] destination = new Double[length];
            System.arraycopy(values, start, destination, offset, count);

            if (selected != null) {
                boolean[] destinationSelected = new boolean[length];
                System.arraycopy(selected, start, destinationSelected, offset, count);
                return getSeletedValues(destination, destinationSelected);
            } else {
                return destination;
            }
        }
    }

    //******************** Trace Integer ****************************
    public class I extends Trace<Integer> {

        public I(String name, int initialSize) {
            super.name = name;
            super.values = new Integer[initialSize];
//            values[0] = initValue; // make getTraceType() working
        }

        public I(String name, Integer[] values) {
            this(name, values.length);
            valueCount = values.length;
            System.arraycopy(values, 0, this.values, 0, values.length);
        }

        public Integer[] getValues(int length, int start, int offset, boolean[] selected) {
            return this.getValues(length, start, offset, valueCount - start, selected);
        }

        public Integer[] getValues(int length, int start, int offset, int count, boolean[] selected) {
            Integer[] destination = new Integer[length];
            System.arraycopy(values, start, destination, offset, count);

            if (selected != null) {
                boolean[] destinationSelected = new boolean[length];
                System.arraycopy(selected, start, destinationSelected, offset, count);
                return getSeletedValues(destination, destinationSelected);
            } else {
                return destination;
            }
        }
    }

    //******************** Trace String ****************************
    public class S extends Trace<String> {

        public S(String name, int initialSize, String initValue) {
            super.name = name;
            if (initialSize > 0) {
                super.values = new String[initialSize];
            }
            values[0] = initValue; // make getTraceType() working
        }

        public S(String name, String[] values) {
            super.name = name;
            super.values = new String[values.length];
            valueCount = values.length;
            System.arraycopy(values, 0, this.values, 0, values.length);
        }

        public String[] getValues(int length, int start, int offset, boolean[] selected) {
            return this.getValues(length, start, offset, valueCount - start, selected);
        }

        public String[] getValues(int length, int start, int offset, int count, boolean[] selected) {
            String[] destination = new String[length];
            System.arraycopy(values, start, destination, offset, count);

            if (selected != null) {
                boolean[] destinationSelected = new boolean[length];
                System.arraycopy(selected, start, destinationSelected, offset, count);
                return getSeletedValues(destination, destinationSelected);
            } else {
                return destination;
            }
        }
    } */

}