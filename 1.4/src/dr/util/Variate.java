/*
 * Variate.java
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

package dr.util;

import dr.stats.DiscreteStatistics;

/**
 * interface for a continuous variate.
 *
 * @version $Id: Variate.java,v 1.15 2006/02/20 17:36:23 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public interface Variate
{
	/**
	 * @return the number of values in this Variate
	 */
	int getCount();

	/**
	 * @return the ith value. Returns NaN if it is missing
	 */
	double get(int index);
	
	/**
	 * @return if the value for index is missing
	 */
	boolean isMissing(int index);
	
	/**
	 * @return the minimum value.
	 */
	double getMin();
	
	/**
	 * @return the maximum value. 
	 */
	double getMax();
	
	/**
	 * @return the range. 
	 */
	double getRange();
	
	/**
	 * @return the sum. 
	 */
	double getSum();
	
	/**
	 * @return the mean. 
	 */
	double getMean();
	
	/**
	 * compute the q-th quantile for the distribution
	 * (= inverse cdf)
	 * 
	 * @param q quantile (0 <= q <= 1)
	 */
	double getQuantile(double q);

	/**
	 * add a value
	 */
	void add(double value);
	
	/**
	 * add all the values
	 */
	void add(Variate values);
	
	/**
	 * add all the values
	 */
	void add(double[] values);
	
	/**
	 * remove a value
	 */
	void remove(int index);
	
	/**
	 * remove all values
	 */
	void removeAll();
	
	/**
	 * A double precision concrete implementation of Variate
	 */
	public class Double implements Variate {
		private double[] values = null;
		private int valueCount = 0;
		
		public Double() {}
		
		public Double(double[] values) {
			this.values = values;
			valueCount = values.length;
		}
		
		/**
		 * add a value
		 */
		public void add(double value) {
		
			if (values == null || valueCount == values.length) {
				double[] newValues = new double[valueCount + 100];
				if (values != null)
					System.arraycopy(values, 0, newValues, 0, values.length);
				values = newValues;
			}
			
			values[valueCount] = value;
			valueCount++;
		}
		
		/**
		 * add all the values
		 */
		public void add(Variate values) {
			for (int i = 0; i < values.getCount(); i++) {
				add(values.get(i));
			}
		}
			
		/**
		 * add all the values
		 */
		public void add(double[] values) {
			for (int i = 0; i < values.length; i++) {
				add(values[i]);
			}
		}
		
	
		/**
		 * remove a value
		 */
		public void remove(int index) {
			if (valueCount > 0) {
				for (int i = index; i < values.length - 1; i++) {
					values[i] = values[i+1];
				}
				valueCount--;
			}
		}
		
		/**
		 * remove all values
		 */
		public void removeAll() { valueCount = 0; }
		
		/**
		 * @return the number of values in this Variate
		 */
		public int getCount() { return valueCount; }

		/**
		 * @return the ith value. 
		 */
		public double get(int index) { return values[index]; }
		
		/**
		 * @return if the value for index is missing
		 */
		public boolean isMissing(int index) { return (values[index] == java.lang.Double.NaN); }
		
		/**
		 * @return the minimum value.
		 */
		public double getMin() {
			double minValue = java.lang.Double.POSITIVE_INFINITY;
			
			for (int i = 0; i < valueCount; i++) {
				if (values[i] < minValue)
					minValue = values[i];
			}
			
			return minValue;
		}
		
		/**
		 * @return the maximum value. 
		 */
		public double getMax() {
			double maxValue = java.lang.Double.NEGATIVE_INFINITY;
			
			for (int i = 0; i < valueCount; i++) {
				if (values[i] > maxValue)
					maxValue = values[i];
			}
			
			return maxValue;
		}
	
		/**
		 * @return the range. 
		 */
		public double getRange() {
			return getMin() - getMax();
		}
		
		/**
		 * @return the mean. 
		 */
		public double getMean() {
			return getSum() / getCount();
		}
		
		/**
		 * @return the sum. 
		 */
		public double getSum() {
			double sum = 0.0;
			
			for (int i = 0; i < getCount(); i++) {
				sum += get(i);
			}
			
			return sum;
		}
		
		/**
		 * compute the q-th quantile for the distribution
		 * (= inverse cdf)
		 * 
		 * @param q quantile (0 <= q <= 1)
		 */
		public double getQuantile(double q) {
			return DiscreteStatistics.quantile(q, values, valueCount);
		}

	};
	
	/**
	 * A single precision concrete implementation of Variate
	 */
	public class Float implements Variate {
		private float[] values = null;
		private int valueCount = 0;
		
		/**
		 * add a value
		 */
		public void add(double value) {
		
			if (values == null || valueCount == values.length) {
				float[] newValues = new float[valueCount + 100];
				if (values != null)
					System.arraycopy(values, 0, newValues, 0, values.length);
				values = newValues;
			}
			
			values[valueCount] = (float)value;
			valueCount++;

		}
		
		/**
		 * add all the values
		 */
		public void add(Variate values) {
			for (int i = 0; i < values.getCount(); i++) {
				add(values.get(i));
			}
		}
			
		/**
		 * add all the values
		 */
		public void add(double[] values) {
			for (int i = 0; i < values.length; i++) {
				add(values[i]);
			}
		}
		
	
		/**
		 * remove a value
		 */
		public void remove(int index) {
			if (valueCount > 0) {
				for (int i = index; i < values.length - 1; i++) {
					values[i] = values[i+1];
				}
				valueCount--;
			}
		}
		
		/**
		 * remove all values
		 */
		public void removeAll() { valueCount = 0; }
		
		
		/**
		 * @return the number of values in this Variate
		 */
		public int getCount() { return valueCount; }

		/**
		 * @return the ith value. 
		 */
		public double get(int index) { return (double)values[index]; }
		
		/**
		 * @return if the value for index is missing
		 */
		public boolean isMissing(int index) { return (values[index] == java.lang.Float.NaN); }
		
		/**
		 * @return the minimum value.
		 */
		public double getMin() {
			double minValue = java.lang.Double.MAX_VALUE;
			
			for (int i = 0; i < valueCount; i++) {
				if (values[i] < minValue)
					minValue = (double)values[i];
			}
			
			return minValue;
		}
		
		/**
		 * @return the maximum value. 
		 */
		public double getMax() {
			double maxValue = java.lang.Double.MIN_VALUE;
			
			for (int i = 0; i < valueCount; i++) {
				if (values[i] > maxValue)
					maxValue = (double)values[i];
			}
			
			return maxValue;
		}
	
		/**
		 * @return the range. 
		 */
		public double getRange() {
			return getMin() - getMax();
		}
		
		/**
		 * @return the mean. 
		 */
		public double getMean() {
			return getSum() / getCount();
		}
		
		/**
		 * @return the sum. 
		 */
		public double getSum() {
			double sum = 0.0;
			
			for (int i = 0; i < valueCount; i++) {
				sum += get(i);
			}
			
			return sum;
		}
		
		/**
		 * compute the q-th quantile for the distribution
		 * (= inverse cdf)
		 * 
		 * @param q quantile (0 <= q <= 1)
		 */
		public double getQuantile(double q) {
			double[] dv = new double[values.length];
			for (int i = 0; i < values.length; i++) {
				dv[i] = (double)values[i];
			}
			return DiscreteStatistics.quantile(q, dv);
		}

	};
}
