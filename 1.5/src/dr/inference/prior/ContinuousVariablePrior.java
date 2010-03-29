/*
 * ContinuousVariablePrior.java
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

//package dr.inference.prior;

/**
 * This class can model flat and scale-independent priors with minimum and maximum hard cutoffs.
 *
 * @author Alexei Drummond
 *
 * @version $Id: ContinuousVariablePrior.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
//public class ContinuousVariablePrior {
//
//	public static final int UNIFORM = 0;
//	public static final int JEFFREYS = 1;
//
//	/**
//	 * Create an uninformative flat prior.
//	 */
//	public ContinuousVariablePrior() {
//	}
//
//	public static ContinuousVariablePrior createJeffreysPrior() {
//		return new ContinuousVariablePrior(0.0, Double.MAX_VALUE, JEFFREYS);
//	}
//
//	public ContinuousVariablePrior(double value) {
//		priorType = UNIFORM;
//		minimumValue = value;
//		maximumValue = value;
//	}
//
//	public ContinuousVariablePrior(double min, double max, int priorType) {
//		minimumValue = min;
//		maximumValue = max;
//		this.priorType = priorType;
//	}
//
//	public void setMinimum( double min_ ){
//		minimumValue = min_;
//	}
//
//	public void setMaximum( double max_ ){
//		maximumValue = max_;
//	}
//
//	public void setPriorType(int priorType) {
//		this.priorType = priorType;
//	}
//
//	public double getMinimum(){
//		return minimumValue;
//	}
//
//	public double getMaximum(){
//		return maximumValue;
//	}
//
//	public int getPriorType() {
//		return priorType;
//	}
//
//	public boolean failed(double value) {
//		return (getLogPrior(value) == Double.NEGATIVE_INFINITY);
//	}
//
//	public double getLogPrior(double value) {
//
//		if ((value < minimumValue) || (value > maximumValue)) {
//
//			return Double.NEGATIVE_INFINITY;
//		}
//
//		if (priorType == JEFFREYS) {
//			return -Math.log(value);
//		} else return 0.0;
//	}
//
//	public boolean fixed() {
//		return minimumValue == maximumValue;
//	}
//
//	public String toString() {
//
//		if (fixed()) {
//			return "value=\"" + minimumValue + "\"";
//		}
//		return "min=\"" + minimumValue + "\" max=\"" + maximumValue + "\" type=\"" + ((priorType == UNIFORM) ? "uniform" : "Jeffreys'") + "\"";
//	}
//
//	private int priorType = UNIFORM;
//	private double minimumValue = -Double.MAX_VALUE;
//	private double maximumValue = Double.MAX_VALUE;
//}
