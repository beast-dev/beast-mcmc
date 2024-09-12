/*
 * MFWithGradient.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.math;


/**
 * interface for a function of several variables with a gradient
 *
 *
 * @author Korbinian Strimmer
 */
public interface MFWithGradient extends MultivariateFunction
{
	/**
	 * compute both function value and gradient at a point
	 *
	 * @param argument  function argument (vector)
	 * @param gradient  gradient (on return)
	 *
	 * @return function value
	 */
	double evaluate(double[] argument, double[] gradient);
	
	/**
	 * compute gradient at a point
	 *
	 * @param argument  function argument (vector)
	 * @param gradient  gradient (on return)
	 */
	void computeGradient(double[] argument, double[] gradient);
}
