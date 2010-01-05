/*
 * UniformParameterPrior.java
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

package dr.inference.prior;

/**
 * This prior provides an uniform and bounded expectation on a parameter.
 */
//public class UniformParameterPrior extends AbstractParameterPrior {
//
//	private double lower;
//	private double upper;
//
//	public UniformParameterPrior(Parameter parameter, double lower, double upper) {
//
//		this(parameter, -1, lower, upper);
//	}
//
//	public UniformParameterPrior(Parameter parameter, int dimension, double lower, double upper) {
//		this.upper = upper;
//		this.lower = lower;
//		setParameter(parameter);
//		setDimension(dimension);
//	}
//
//	public final double getLogPriorComponent(double value) {
//		if (value >= lower && value <= upper) return 0.0;
//		return Double.NEGATIVE_INFINITY;
//	}
//
//	public Element createElement(Document d) {
//		Element e = d.createElement("uniformPrior");
//		e.setAttribute("lower", lower + "");
//		e.setAttribute("upper", upper + "");
//		return e;
//	}
//
//	public double getLowerLimit() { return lower; }
//	public double getUpperLimit() { return upper; }
//
//	public String toString() {
//
//		StringBuilder buffer = new StringBuilder();
//		if (lower == -Double.MAX_VALUE) {
//            buffer.append("(").append(formatter.format(Double.NEGATIVE_INFINITY).trim());
//		} else if (lower == Double.MIN_VALUE) {
//            buffer.append("(").append(formatter.format(0.0).trim());
//		} else {
//            buffer.append("[").append(formatter.format(lower).trim());
//		}
//
//		buffer.append(", ");
//
//		if (upper == Double.MAX_VALUE) {
//            buffer.append(formatter.format(Double.POSITIVE_INFINITY).trim()).append(")");
//		} else if (upper == -Double.MIN_VALUE) {
//            buffer.append(formatter.format(0.0).trim()).append(")");
//		} else {
//            buffer.append(formatter.format(upper).trim()).append("]");
//		}
//
//		return buffer.toString();
//	}
//
//	public String toHTML() {
//		return "<font color=\"#FF00FF\">" + toString() + "</font>";
//	}
//}
