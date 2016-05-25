/*
 * MultivariateIntegral.java
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

package dr.math;

/**
 * An interface for integration of univariate functions.
 *
 * @author Alexei Drummond
 *
 * @version $Id: Integral.java,v 1.3 2005/05/24 20:26:01 rambaut Exp $
 */
public interface MultivariateIntegral {

    /**
     * @return the approximate integral of the given function
     * within the given range.
     * @param f the function whose integral is of interest
     * @param mins the minimum values of the function in each dimension
     * @param maxes the upper limit of the function in each dimension
     */
    public double integrate(MultivariateFunction f, double[] mins, double[] maxes);
}
