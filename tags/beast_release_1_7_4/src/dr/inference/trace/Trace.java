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
import java.util.Collections;
import java.util.List;

/**
 * A simple class that stores a trace for a single statistic
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Trace.java,v 1.11 2005/07/11 14:07:26 rambaut Exp $
 */
public class Trace<T> {

//    public static final int INITIAL_SIZE = 1000;
//    public static final int INCREMENT_SIZE = 1000;

    // use <Double> for integer, but traceType must = INTEGER
    protected TraceFactory.TraceType traceType = TraceFactory.TraceType.DOUBLE;
    protected List<T> values = new ArrayList<T>();
    //    protected int valueCount = 0;
    protected String name;

//    private Object[] range;

    public Trace(String name) { // traceType = TraceFactory.TraceType.DOUBLE; 
        this.name = name;
    }

    public Trace(String name, TraceFactory.TraceType traceType) {
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

    public int getValuesSize() {
        return values.size();
    }

    public T getValue(int index) {
        return values.get(index); // filter?
    }

    public String[] getRange() { // Double => bounds; Integer and String => unique values
        String[] range;
        if (getTraceType() == TraceFactory.TraceType.DOUBLE) {
            range = new String[2];
            range[0] = Double.toString(Double.MAX_VALUE);
            range[1] = Double.toString(Double.MIN_VALUE);
            for (Object t : values) {
                if (Double.parseDouble(range[0]) < (Double) t) range[0] = t.toString();
                if (Double.parseDouble(range[1]) > (Double) t) range[1] = t.toString();
            }

        } else {
            List<String> r = new ArrayList<String>();
            for (Object t : values) {
                if (traceType == TraceFactory.TraceType.INTEGER) { // as Integer is stored as Double in Trace
                    if (!r.contains(Integer.toString(((Number) t).intValue())))
                        r.add(Integer.toString(((Number) t).intValue()));
                } else {
                    if (!r.contains(t.toString()))
                        r.add(t.toString());
                }
            }

            Collections.sort(r);
            range = new String[r.size()];
            range = r.toArray(range);

        }

        return range;
    }

    /**
     * @param fromIndex low endpoint (inclusive) of the subList.
     * @param toIndex   high endpoint (exclusive) of the subList.
     * @param selected  if null then no filter, otherwise selected.length should = getValuesSize()
     * @return The list of values (which are selected values if filter applied)
     */
    public List<T> getValues(int fromIndex, int toIndex, boolean[] selected) {
        if (toIndex > getValuesSize() || fromIndex > toIndex)
            throw new RuntimeException("Invalid index : fromIndex = " + fromIndex + "; toIndex = " + toIndex
                    + "; List size = " + getValuesSize() + "; in Trace " + name);

        if (selected != null) {
            if (selected.length != getValuesSize())
                throw new RuntimeException("size of values is different with selected[] in Trace " + name);

            List<T> valuesList = new ArrayList<T>();
            for (int i = fromIndex; i < toIndex; i++) {
                if (selected[i])
                    valuesList.add(values.get(i));
            }
            if (valuesList.size() < 1) throw new RuntimeException("There is no value sent by getValue() !");
            return valuesList;
        } else {
            return values.subList(fromIndex, toIndex);
        }
    }

    public String getName() {
        return name;
    }

//    public Class getTraceType() {
//        if (values.get(0) == null) {
//            return null;
//        }
//        return values.get(0).getClass();
//    }

    public TraceFactory.TraceType getTraceType() {
        return traceType;
    }

    public void setTraceType(TraceFactory.TraceType traceType) {
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

    public boolean isIn(int i) {
        return filter.isIn(values.get(i), traceType);
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