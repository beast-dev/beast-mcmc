/*
 * AbstractParameterPrior.java
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
 * This class provides an abstract superclass for all parameter priors.
 */
//public abstract class AbstractParameterPrior implements ParameterPrior {
//
//	//************************************************************************
//	// abstract methods
//	//************************************************************************
//
//	/**
//	 * @return the log prior of a single-value component of the parameter.
//	 */
//	public abstract double getLogPriorComponent(double value);
//
//	//************************************************************************
//	// final methods
//	//************************************************************************
//
//	public final void setParameter(Parameter param) { this.parameter = param; }
//	public final Parameter getParameter() { return parameter; }
//	public final void setDimension(int dim) { dimension = dim; }
//
//	public final double getLogPrior(Model model) {
//
//		if (dimension == -1) {
//			double logL = 0.0;
//			for (int i =0; i < parameter.getDimension(); i++) {
//				logL += getLogPriorComponent(parameter.getParameterValue(i));
//			}
//			return logL;
//		} else return getLogPriorComponent(parameter.getParameterValue(dimension));
//	}
//
//	public final String getPriorName() { return toString(); }
//
//	//************************************************************************
//	// protected instance variables
//	//************************************************************************
//
//	protected NumberFormatter formatter = new NumberFormatter(6);
//
//	//************************************************************************
//	// private instance variables
//	//************************************************************************
//
//	/** the parameter this prior acts on */
//	private Parameter parameter = null;
//
//	/** the dimension of the parameter that this prior works on, -1 signifies all dimensions */
//	private int dimension = -1;
//}
