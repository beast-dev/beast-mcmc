/*
 * Parameter.java
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

package dr.inference.model;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;

/**
 * Represents a multi-dimensional continuous parameter. 
 *
 * @version $Id: Parameter.java,v 1.22 2005/06/08 11:23:25 alexei Exp $
 *
 * @author Alexei Drummond
 */
public interface Parameter extends Statistic
{
	/** @return the parameter's scalar value in the given dimension
     *  @param dim
     */
	double getParameterValue(int dim);

    /** @return the parameter's values
     */
    double[] getParameterValues();

	/** sets the scalar value in the given dimension of this parameter
     * @param dim
     * @param value
     */
	void setParameterValue(int dim, double value);

    /** sets the scalar value in the given dimensin of this parameter to val, without firing any events
     * @param dim
     * @param value
     */
    void setParameterValueQuietly(int dim, double value);

	/** @return the name of this parameter */
	String getParameterName();

	/** adds a parameter listener that is notified when this parameter changes.
     * @param listener
     */
	void addParameterListener(ParameterListener listener);
	
	/** removes a parameter listener.
     * @param listener\
     */
	void removeParameterListener(ParameterListener listener);
	
	/** stores the state of this parameter for subsquent restore */
	void storeParameterValues();
	
	/** restores the stored state of this parameter */
	void restoreParameterValues();
	
	/** accepts the stored state of this parameter */
	void acceptParameterValues();
	
	/** adopt the state of the source parameter
     * @param source
     */
	void adoptParameterValues(Parameter source);
	
	/**
     * @return true if values in all dimensions are within their bounds
     * */
	boolean isWithinBounds();
	
	/** 
     * Can be called before store is called. If it results in new
	 * dimensions, then the value of the first dimension is copied into the new dimensions.
     * @param dim new dimention
     */
	public void setDimension(int dim);
	
	/** Adds new bounds to this parameter
     * @param bounds to add
     */
	void addBounds(Bounds bounds);
	
	/** @return the intersection of all bounds added to this parameter */
	Bounds getBounds();
		
	/**
	 * Abstract base class for parameters
	 */
	public abstract class Abstract extends Statistic.Abstract implements Parameter {
		
		public int getDimension() { return 1; }
		
		/**
		 * Fired when all dimensions of the parameter have changed
		 */
		public void fireParameterChangedEvent() {
			fireParameterChangedEvent(-1);
		}
		
		/**
		 * Fired when a single dimension of the parameter has changed
         * @param index which dimention changed
         */
		public void fireParameterChangedEvent(int index) {
			if (listeners != null) {
                for (ParameterListener listener : listeners) {
                    listener.parameterChangedEvent(this, index);
                }
            }
		}
		
		public void addParameterListener(ParameterListener listener) {
			if (listeners == null) {listeners = new ArrayList<ParameterListener>();}
			listeners.add(listener);
		}
		
		public void removeParameterListener(ParameterListener listener) {
			if (listeners != null) {
				listeners.remove(listener);
			}
		}
		
		public final String getStatisticName() { return getParameterName(); }
		public final double getStatisticValue(int dim) { return getParameterValue(dim); }
		
		public void setDimension(int dim) {
			throw new UnsupportedOperationException();
		}

        /**
         * Defensively returns copy of parameter array.
         * @return
         */
        public double[] getParameterValues() {

            double[] copyOfValues = new double[getDimension()];
            for (int i = 0; i < copyOfValues.length; i++) {
                copyOfValues[i] = getParameterValue(i);
            }
            return copyOfValues;
        }


		public final void storeParameterValues() {
			if (isValid) {
				storeValues();
				
				isValid = false;
			}
		}
		
		public final void restoreParameterValues() {
			if (!isValid) {
				restoreValues();
				
				isValid = true;
			}
		}
		
		public final void acceptParameterValues() {
			if (!isValid) {
				acceptValues();

				isValid = true;
			}
		}
		
		public final void adoptParameterValues(Parameter source) {
		
			adoptValues(source);

			isValid = true;
		}
		
		public boolean isWithinBounds() {
			Bounds bounds = getBounds();
			for (int i = 0; i < getDimension(); i++) {
				if (getParameterValue(i) < bounds.getLowerLimit(i) ||
					getParameterValue(i) > bounds.getUpperLimit(i)) {
					return false;
				}
			}
			return true;
		}
	
		protected abstract void storeValues();
		protected abstract void restoreValues();
		protected abstract void acceptValues();
		protected abstract void adoptValues(Parameter source);

		public String toString() {
			StringBuffer buffer = new StringBuffer(String.valueOf(getParameterValue(0)));
            buffer.append(getId()).append("=[").append(String.valueOf(getBounds().getLowerLimit(0)));
            buffer.append(",").append(String.valueOf(getBounds().getUpperLimit(0))).append("]");
			
			for (int i = 1; i < getDimension(); i++) {
                buffer.append(", ").append(String.valueOf(getParameterValue(i)));
                buffer.append("[").append(String.valueOf(getBounds().getLowerLimit(i)));
                buffer.append(",").append(String.valueOf(getBounds().getUpperLimit(i))).append("]");
			}
			return buffer.toString();
		}
		
		public Element createElement(Document document) {
			throw new IllegalArgumentException();
		}
		
		private boolean isValid = true;
		
		private ArrayList<ParameterListener> listeners;
	}
	
	
	/**
	 * A class that implements the Parameter interface.
	 */
	class Default extends Abstract {
	
		public Default(int dimension) {
			this(dimension, 1.0);
		}
	
		public Default(double initialValue) {
			values = new double[1];
			values[0] = initialValue;
			this.bounds = null;
		}
		
		public Default(int dimension, double initialValue) {
			values = new double[dimension];
			for (int i =0; i < dimension; i++) {
				values[i] = initialValue;
			} 
			this.bounds = null;
		}
		
		public Default(double[] values) {
			this.values = new double[values.length];
			for (int i =0; i < values.length; i++) {
				this.values[i] = values[i];
			}
		}
		
		public void addBounds(Bounds boundary) {
			if (bounds == null) {
				bounds = new IntersectionBounds(getDimension());
			}
			bounds.addBounds(boundary);
			
			// can't change dimension after bounds are added!
			hasBeenStored = true;
		}
			
		//********************************************************************
		// GETTERS
		//********************************************************************
		
		public int getDimension() { return values.length; }
		
		public double getParameterValue(int i) { return values[i]; }

        /**
         * Defensively returns copy of parameter array.
         * @return
         */
        public final double[] getParameterValues() {

            double[] copyOfValues = new double[values.length];
            System.arraycopy(values,0,copyOfValues,0,values.length);
            return copyOfValues;
        }

		public Bounds getBounds() { 
			if (bounds == null) {
//				bounds = new IntersectionBounds(getDimension());
				throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
			}
			return bounds; 
		}
		
		public String getParameterName() { return getId(); }
		
		//********************************************************************
		// SETTERS
		//********************************************************************
		
		/** 
		 * Can only be called before store is called. If it results in new
		 * dimensions, then the value of the first dimension is copied into the new dimensions.
		 */
		public void setDimension(int dim) {
			if (!hasBeenStored) {
				double[] newValues = new double[dim];
				for (int i = 0; i < values.length; i++) {
					newValues[i] = values[i];
				}
				for (int i = values.length; i < newValues.length; i++) {
					newValues[i] = values[0];
				}
				values = newValues;
			} else throw new RuntimeException("Can't change dimension after store has been called!");
		}
		
		public void setParameterValue(int i, double val) { 
			values[i] = val;
			
			fireParameterChangedEvent(i);
		}

        /**
         * Sets the value of the parameter without firing a changed event.
         * @param dim
         * @param value
         */
        public void setParameterValueQuietly(int dim, double value) {
            values[dim] = value;
        }

		protected final void storeValues() {
			hasBeenStored = true;
			if (storedValues == null) {
				storedValues = new double[values.length];
			}
			System.arraycopy(values, 0, storedValues, 0, values.length);
		}
		
		protected final void restoreValues() {
			
			//swap the arrays
			double[] temp = storedValues;
			storedValues = values;
			values = temp;
			
			//if (storedValues != null) {
			//	System.arraycopy(storedValues, 0, values, 0, values.length);
			//} else throw new RuntimeException("restore called before store!");
		}
		
		/** Nothing to do */
		protected final void acceptValues() { }
			
		protected final void adoptValues(Parameter source) {
		
			if (getDimension() != source.getDimension()) {
				throw new RuntimeException("The two parameters don't have the same number of dimensions");
			}
			
			for (int i = 0, n = getDimension(); i < n; i++) {
				values[i] = source.getParameterValue(i);
			}
		}
		
		private double[] values;
		
		private double[] storedValues;
		
		private boolean hasBeenStored = false;
		private IntersectionBounds bounds = null;
	}

	class DefaultBounds implements Bounds {
		
		public DefaultBounds(double upper, double lower, int dimension) {
			
			this.uppers = new double[dimension];
			this.lowers = new double[dimension];
			for (int i =0; i < dimension; i++) {
				uppers[i] = upper;
				lowers[i] = lower;
			}
		}
		
		public DefaultBounds(ArrayList<Double> upperList, ArrayList<Double> lowerList) {

            final int length = upperList.size();
            if (length != lowerList.size()) {
				throw new IllegalArgumentException("upper and lower limits must be defined on the same number of dimensions.");
			}
			uppers = new double[length];
			lowers = new double[length];
			for (int i = 0; i < uppers.length; i++) {
				uppers[i] = upperList.get(i);
				lowers[i] = lowerList.get(i);
			}
		}
		
		public DefaultBounds(double[] uppers, double[] lowers) {
			
			if (uppers.length != lowers.length) { 
				throw new IllegalArgumentException("upper and lower limits must be defined on the same number of dimensions.");
			}
			this.uppers = uppers;
			this.lowers = lowers;
		}
		
		public double getUpperLimit(int i) { return uppers[i]; }
		public double getLowerLimit(int i) { return lowers[i]; }
		public int getBoundsDimension() { return uppers.length; }	
	
		private double[] uppers, lowers;
	}
}
