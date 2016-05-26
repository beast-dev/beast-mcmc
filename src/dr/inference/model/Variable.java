/*
 * Variable.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.util.Identifiable;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic random variable.
 *
 * @author Alexei Drummond
 */
public interface Variable<V> extends Identifiable {

    public enum ChangeType {
        VALUE_CHANGED,
        REMOVED,
        ADDED,
        ALL_VALUES_CHANGED
    }

    /**
     * @return the name of this variable.
     */
    public String getVariableName();

    public V getValue(int index);

    public void setValue(int index, V value);

    public V[] getValues();

    /**
     * @return the size of this variable - i.e. the length of the vector
     */
    public int getSize();

    /**
     * adds a parameter listener that is notified when this parameter changes.
     *
     * @param listener the listener
     */
    void addVariableListener(VariableListener listener);

    /**
     * removes a parameter listener.
     *
     * @param listener the listener
     */
    void removeVariableListener(VariableListener listener);

    /**
     * stores the state of this parameter for subsequent restore
     */
    void storeVariableValues();

    /**
     * restores the stored state of this parameter
     */
    void restoreVariableValues();

    /**
     * accepts the stored state of this parameter
     */
    void acceptVariableValues();

    /**
     * @return the bounds on this parameter
     */
    Bounds<V> getBounds();

    void addBounds(Bounds<V> bounds);

    public abstract class Base<V> implements Variable<V>, Loggable {
        Base(String id) {
            this.id = id;
        }

        protected void fireVariableChanged(int index) {
            for (VariableListener listener : listeners) {
                listener.variableChangedEvent(this, index, ChangeType.VALUE_CHANGED);
            }
        }

        public void addVariableListener(VariableListener listener) {
            listeners.add(listener);
        }

        public void removeVariableListener(VariableListener listener) {
            listeners.remove(listener);
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }

        public String getVariableName() {
            return id;
        }

        protected List<VariableListener> listeners = new ArrayList<VariableListener>();

        protected String id;
    }


    public abstract class BaseNumerical<V extends Number> extends Base<V> {
        BaseNumerical(String id) {
            super(id);
        }

        class StatisticColumn extends NumberColumn {
            private final int dim;

            public StatisticColumn(String label, int dim) {
                super(label);
                this.dim = dim;
            }

            public double getDoubleValue() {
                return getValue(dim).doubleValue();
            }
        }

        /**
         * @return the log columns.
         */
        public LogColumn[] getColumns() {
            LogColumn[] columns = new LogColumn[getSize()];
            if (getSize() == 1) {
                columns[0] = new StatisticColumn(this.getVariableName(), 0);
            } else {
                for (int i = 0; i < getSize(); i++) {
                    columns[i] = new StatisticColumn(this.getVariableName() + "[" + i + "]", i);
                }
            }
            return columns;
        }
    }

    public class D implements Variable<Double>, Loggable {

        public D(double value, int size) {
            values = new double[size];
            storedValues = new double[values.length];

            for (int i = 0; i < size; i++) {
                values[i] = value;
            }
        }

        public D(double[] v) {
            values = new double[v.length];
            System.arraycopy(v, 0, values, 0, v.length);
            storedValues = new double[values.length];
        }

        public D(String name, double[] v) {
            this(v);
            setId(name);
        }

        public D(String name, double value) {
            this(name, new double[]{value});
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }


        public String getVariableName() {
            return id;
        }

        public Double getValue(int index) {
            return values[index];
        }

        public Double[] getValues() {
            Double[] copyOfValues = new Double[values.length];
            for (int i = 0; i < values.length; i++) {
                copyOfValues[i] = values[i];
            }
            return copyOfValues;
        }

        public void setValue(int index, Double value) {
            values[index] = value;
            fireVariableChanged(index);
        }

        private void fireVariableChanged(int index) {
            for (VariableListener listener : listeners) {
                listener.variableChangedEvent(this, index, ChangeType.VALUE_CHANGED);
            }
        }

        public int getSize() {
            return values.length;
        }

        public void addVariableListener(VariableListener listener) {
            listeners.add(listener);
        }

        public void removeVariableListener(VariableListener listener) {
            listeners.remove(listener);
        }

        public void storeVariableValues() {
            System.arraycopy(values, 0, storedValues, 0, storedValues.length);
        }

        public void restoreVariableValues() {
            double[] temp = storedValues;
            storedValues = values;
            values = temp;
        }

        public void acceptVariableValues() {
        }

        public Bounds<Double> getBounds() {
            if (bounds == null) {
                return new Bounds<Double>() {
                    public Double getUpperLimit(int dimension) {
                        return Double.MAX_VALUE;
                    }

                    public Double getLowerLimit(int dimension) {
                        return -Double.MAX_VALUE;
                    }

                    public int getBoundsDimension() {
                        return getSize();
                    }
                };
            } else return bounds;
        }

        // **************************************************************
        // Loggable IMPLEMENTATION
        // **************************************************************

        /**
         * @return the log columns.
         */
        public LogColumn[] getColumns() {
            LogColumn[] columns = new LogColumn[getSize()];
            if (getSize() == 1) {
                columns[0] = new StatisticColumn(getVariableName(), 0);
            } else {
                for (int i = 0; i < getSize(); i++) {
                    columns[i] = new StatisticColumn(getVariableName() + "[" + i + "]", i);
                }
            }
            return columns;
        }

        /**
         * Careful use! Do not write to the array provided!!
         *
         * @return the underlying array of doubles
         */
        public double[] peekValues() {
            return values;
        }

        public void addBounds(Bounds<Double> b) {
            if (bounds == null) {
                bounds = new IntersectionBounds(getSize());
            }
            bounds.addBounds(b);
        }

        private class StatisticColumn extends NumberColumn {
            private final int dim;

            public StatisticColumn(String label, int dim) {
                super(label);
                this.dim = dim;
            }

            public double getDoubleValue() {
                return getValue(dim);
            }
        }

        String id;
        double[] values;
        double[] storedValues;
        List<VariableListener> listeners = new ArrayList<VariableListener>();
        private IntersectionBounds bounds = null;
    }

    public class DM implements Variable<double[]>, Loggable {

        public DM(double[][] v) {

            lower = new double[v.length];
            upper = new double[v.length];

            values = new double[v.length][v[0].length];
            for (int i = 0; i < v.length; i++) {
                System.arraycopy(v[i], 0, values[i], 0, v[i].length);
                lower[i] = -Double.MAX_VALUE;
                upper[i] = Double.MAX_VALUE;
            }
            storedValues = new double[values.length][values[0].length];
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }


        public String getVariableName() {
            return id;
        }

        public double[] getValue(int index) {
            return values[index];
        }

        public double[][] getValues() {
            double[][] copyOfValues = new double[values.length][values[0].length];
            for (int i = 0; i < values.length; i++) {
                System.arraycopy(values[i], 0, copyOfValues[i], 0, values[i].length);
            }
            return copyOfValues;
        }

        public void setValue(int index, double[] value) {
            System.arraycopy(value, 0, values[index], 0, values[index].length);

            fireVariableChanged(index);
        }

        private void fireVariableChanged(int index) {
            for (VariableListener listener : listeners) {
                listener.variableChangedEvent(this, index, ChangeType.VALUE_CHANGED);
            }
        }

        public int getSize() {
            return values.length;
        }

        public void addVariableListener(VariableListener listener) {
            listeners.add(listener);
        }

        public void removeVariableListener(VariableListener listener) {
            listeners.remove(listener);
        }

        public void storeVariableValues() {
            for (int i = 0; i < values.length; i++) {
                System.arraycopy(values[i], 0, storedValues[i], 0, storedValues[i].length);
            }
        }

        public void restoreVariableValues() {
            double[][] temp = storedValues;
            storedValues = values;
            values = temp;
        }

        public void acceptVariableValues() {
        }

        public Bounds<double[]> getBounds() {
            return new Bounds<double[]>() {
                public double[] getUpperLimit(int dimension) {
                    return upper;
                }

                public double[] getLowerLimit(int dimension) {
                    return lower;
                }

                public int getBoundsDimension() {
                    return getSize();
                }
            };
        }

        public void addBounds(Bounds<double[]> bounds) {

        }


        // **************************************************************
        // Loggable IMPLEMENTATION
        // **************************************************************

        /**
         * @return the log columns.
         */
        public LogColumn[] getColumns() {
            LogColumn[] columns = new LogColumn[getSize()];
            for (int i = 0; i < getSize(); i++) {
                double[] values = getValue(i);
                for (int j = 0; j < values.length; j++) {
                    if (getSize() == 1) {
                        columns[i] = new StatisticColumn(getVariableName() + "[" + j + "]", i, j);
                    } else {
                        columns[i] = new StatisticColumn(getVariableName() + "[" + i + "," + j + "]", i, j);
                    }
                }
            }
            return columns;
        }

        private class StatisticColumn extends NumberColumn {
            private final int i;
            private final int j;

            public StatisticColumn(String label, int i, int j) {
                super(label);
                this.i = i;
                this.j = j;
            }

            public double getDoubleValue() {
                return getValue(i)[j];
            }
        }

        String id;
        double[][] values;
        double[][] storedValues;
        List<VariableListener> listeners = new ArrayList<VariableListener>();
        double[] lower;
        double[] upper;
    }

    public class I implements Variable<Integer>, Loggable {

        public I(int value, int size) {
            values = new int[size];
            storedValues = new int[values.length];

            for (int i = 0; i < size; i++) {
                values[i] = value;
            }
        }

        public I(int[] v) {
            values = new int[v.length];
            System.arraycopy(v, 0, values, 0, v.length);
            storedValues = new int[values.length];
        }

        public I(String name, int[] v) {
            this(v);
            setId(name);
        }

        public I(String name, int value) {
            this(name, new int[]{value});
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }


        public String getVariableName() {
            return id;
        }

        public Integer getValue(int index) {
            return values[index];
        }

        public Integer[] getValues() {
            Integer[] copyOfValues = new Integer[values.length];
            for (int i = 0; i < values.length; i++) {
                copyOfValues[i] = values[i];
            }
            return copyOfValues;
        }

        public void setValue(int index, Integer value) {
            values[index] = value;
            fireVariableChanged(index);
        }

        private void fireVariableChanged(int index) {
            for (VariableListener listener : listeners) {
                listener.variableChangedEvent(this, index, ChangeType.VALUE_CHANGED);
            }
        }

        public int getSize() {
            return values.length;
        }

        public void addVariableListener(VariableListener listener) {
            listeners.add(listener);
        }

        public void removeVariableListener(VariableListener listener) {
            listeners.remove(listener);
        }

        public void storeVariableValues() {
            System.arraycopy(values, 0, storedValues, 0, storedValues.length);
        }

        public void restoreVariableValues() {
            int[] temp = storedValues;
            storedValues = values;
            values = temp;
        }

        public void acceptVariableValues() {
        }

        public Bounds<Integer> getBounds() {
            if (bounds == null) {
                return new Bounds<Integer>() {
                    public Integer getUpperLimit(int dimension) {
                        return Integer.MAX_VALUE;
                    }

                    public Integer getLowerLimit(int dimension) {
                        return -Integer.MAX_VALUE;
                    }

                    public int getBoundsDimension() {
                        return getSize();
                    }
                };
            } else return bounds;
        }

        // **************************************************************
        // Loggable IMPLEMENTATION
        // **************************************************************

        /**
         * @return the log columns.
         */
        public LogColumn[] getColumns() {
            LogColumn[] columns = new LogColumn[getSize()];
            if (getSize() == 1) {
                columns[0] = new StatisticColumn(getVariableName(), 0);
            } else {
                for (int i = 0; i < getSize(); i++) {
                    columns[i] = new StatisticColumn(getVariableName() + "[" + i + "]", i);
                }
            }
            return columns;
        }

        /**
         * Careful use! Do not write to the array provided!!
         *
         * @return the underlying array of doubles
         */
        public int[] peekValues() {
            return values;
        }

        public void addBounds(Bounds<Integer> b) {
            if (bounds == null) {
                bounds = b;
            } else {
                if (!(bounds instanceof Bounds.Staircase)) {
                    Bounds.Staircase newBounds = new Bounds.Staircase(getSize());
                    newBounds.addBounds(bounds);
                    bounds = newBounds;
                }

                ((Bounds.Staircase) bounds).addBounds(b);
            }

        }

        private class StatisticColumn extends NumberColumn {
            private final int dim;

            public StatisticColumn(String label, int dim) {
                super(label);
                this.dim = dim;
            }

            public double getDoubleValue() {
                return getValue(dim);
            }
        }

        String id;
        int[] values;
        int[] storedValues;
        List<VariableListener> listeners = new ArrayList<VariableListener>();
        private Bounds<Integer> bounds = null;
    }
}
