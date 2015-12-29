/*
 * VonMisesDistribution.java
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

import cern.jet.math.Bessel;
import dr.geo.math.Space;

/**
 * @author Marc Suchard
 */
public class VonMisesDistribution extends HyperSphereDistribution {
        
    public VonMisesDistribution(int dim, Space space, double[] mean, double kappa) {
        super(dim, space, mean, kappa);
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, kappa, space);
    }

    public String getType() {
        return "von Mises";
    }

    protected int getAllowableDim() {
        return 2;
    }

    private static double logNormalizationConstant(double kappa) {
         return -LOG_2_PI - Math.log(Bessel.i0(kappa));
    }

    public static double logPdf(double[] x, double[] mean, double kappa, Space space) {
        return logNormalizationConstant(kappa) + kappa * HyperSphereDistribution.innerProduct(x, mean, space);
    }

    public static void main(String[] arg) {        
        // Test in radians
        double kappa = 2;
        double[] mean = { 1 };
        double[] x = { -4 };
        System.err.println("logP = "+logPdf(x, mean, kappa, Space.RADIANS)+" ?= -2.094546");
    }
}