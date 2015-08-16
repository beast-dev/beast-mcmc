/*
 * HalfTDistribution.java
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

package dr.math.distributions;

/**
 * @author Marc A. Suchard
 */
public class HalfTDistribution extends TDistribution {

    public HalfTDistribution(double scale, double df) {
        super(0.0, scale, df);
    }

    public double pdf(double x) {
        return x < 0.0 ? 0.0 : super.pdf(x) * 2.0;
    }

    public double logPdf(double x) {
        return x < 0.0 ? Double.NEGATIVE_INFINITY : super.logPdf(x) + Math.log(2.0);
    }

    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

}
