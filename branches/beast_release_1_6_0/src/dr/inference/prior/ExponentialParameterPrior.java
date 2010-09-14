/*
 * ExponentialParameterPrior.java
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
//
//import dr.inference.model.Parameter;
//import org.w3c.dom.Document;
//import org.w3c.dom.Element;

/**
 * This prior provides an exponential expectation on a single parameter
 * @author Alexei Drummond
 * @version $Id: ExponentialParameterPrior.java,v 1.5 2005/05/24 20:26:00 rambaut Exp $
 */
//public class ExponentialParameterPrior extends AbstractParameterPrior {
//
//	/**
//	 * the mean of the exponential distribution.
//	 */
//	double mean;
//
//	public ExponentialParameterPrior(Parameter parameter, double mean) {
//		this(parameter, -1, mean);
//	}
//
//	public ExponentialParameterPrior(Parameter parameter, int dimension, double mean) {
//		this.mean = mean;
//		setParameter(parameter);
//		setDimension(dimension);
//	}
//
//	public double getLogPriorComponent(double value) {
//		return - value / mean;
//	}
//
//	public Element createElement(Document d) {
//		Element e = d.createElement("exponentialPrior");
//		e.setAttribute("mean", mean + "");
//		return e;
//	}
//
//	public final double getMean() { return mean; }
//
//	public String toString() {
//		return "Exponential(" + formatter.format(mean).trim() + ")";
//	}
//
//	public String toHTML() {
//		return "<font color=\"#FF00FF\">Exponential(" + formatter.format(mean).trim() + ")</font>";
//	}
//}
