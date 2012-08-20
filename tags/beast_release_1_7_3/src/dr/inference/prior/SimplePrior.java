/*
 * SimplePrior.java
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

/**
 * This class provides an abstract superclass for simple priors of single variables. 
 * The prior has an upper and lower limit and uniform or Jeffery's prior
 * within these boundaries.
 *
 * @author Alexei Drummond
 *
 * @version $Id: SimplePrior.java,v 1.8 2005/05/24 20:26:00 rambaut Exp $
 */
//public class SimplePrior implements Prior {
//
//	ContinuousVariablePrior cvp;
//	Parameter parameter;
//
//	public SimplePrior(ContinuousVariablePrior cvp, Parameter parameter) {
//		this.cvp = cvp;
//		this.parameter = parameter;
//	}
//
//	public double getLogPrior(dr.inference.model.Model model) {
//
//		double logPrior = 0.0;
//		for (int i =0; i < parameter.getDimension(); i++) {
//			logPrior += cvp.getLogPrior(parameter.getParameterValue(i));
//		}
//
//		return logPrior;
//	}
//
//	public final String getPriorName() {
//		return parameter.getParameterName();
//	}
//}
