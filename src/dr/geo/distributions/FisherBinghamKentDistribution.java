/*
 * FisherBinghamKentDistribution.java
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

package dr.geo.distributions;

import dr.geo.math.Space;

/**
 * @author Marc Suchard
 */
public class FisherBinghamKentDistribution extends HyperSphereDistribution {

    public FisherBinghamKentDistribution(int dim, Space space, double[] mean, double[] major, double[] minor,
                                         double kappa, double beta) {
        super(dim, space, mean, kappa);
        this.major = major;
        this.minor = minor;
        this.beta = beta;
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, major, minor);
    }

    public static double logPdf(double[] x, double[] gamma1, double[] gamma2, double[] gamma3) {
        return 0;
    }

    public String getType() {
        return "Fisher-Bingham (Kent)";
    }

    protected int getAllowableDim() {
        return 3;
    }

    private double[] major;
    private double[] minor;
    private double beta;
}