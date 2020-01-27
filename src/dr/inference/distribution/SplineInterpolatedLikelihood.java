/*
 * SplineInterpolatedLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.SplineBasis;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class SplineInterpolatedLikelihood extends EmpiricalDistributionLikelihood {

    private static final double outsideLogDensity = Double.NEGATIVE_INFINITY; // Use for a proper posterior

    public SplineInterpolatedLikelihood(String fileName, int degree, boolean inverse, boolean byColumn) {
        super(fileName, inverse, byColumn);

         // Set-up spline basis, could be degree = 1 for linear interpolation
//        splineBasis = new SplineBasis(getId(),new Variable.D(values), new Variable.D(density), degree);

        // Something is wrong with the spline basis routines...  just do simple linear interpolation

    }

    public SplineInterpolatedLikelihood(List<EmpiricalDistributionData> dataList, int degree, boolean inverse) {
        super(dataList, inverse);
    }

    @Override
    protected double logPDF(double x, EmpiricalDistributionData data) {
//        return splineBasis.evaluate(x);

        final double[] values = data.values;
        final double[] density = data.density;
        final int len = values.length;

        if (x < values[0] || x > values[len - 1])
            return outsideLogDensity;

        double rtnValue = 0;

        // Binary search for x below beta:
//        unsigned int low = 0;
//        unsigned int high = x.size() - 1;
//        unsigned int m = 0;
//        while (high > low) {
//            m = floor(low + (high-low) / 2.0);
//            if (x[m] > beta[i] || low == m)
//                high = m;
//            else
//                low = m;
//        }

        for(int i=1; i<len; i++) {
            if (values[i] > x) { // first largest point
                final double diffValue = values[i] - values[i-1];
                final double diffDensity = density[i] - density[i-1];
                rtnValue = density[i] - (values[i]-x) / diffValue * diffDensity;
                break;
            }
        }

        if (!data.densityInLogSpace) {
            rtnValue = Math.log(rtnValue);
        }

        if (inverse)
            rtnValue *= -1;
  
        return rtnValue;
    }

    @Override
    protected double gradientLogPdf(double x, EmpiricalDistributionData data) {

        final double[] values = data.values;
        final double[] density = data.density;
        final int len = values.length;

        if (x < values[0] || x > values[len - 1]) {
            return 0.0;
        }

        double gradient = 0;

        for (int i = 0; i < len; ++i) {
            if (values[i] > x) {
                final double diffValue = values[i] - values[i - 1];
                final double diffDensity = density[i] - density[i - 1];
                gradient = diffDensity / diffValue;
                break;
            }
        }

        if (!data.densityInLogSpace) {
            throw new RuntimeException("Not yet implemented"); // TODO
        }

        if (inverse) {
            throw new RuntimeException("Not yet implemented"); // TODO
        }

        return gradient;
    }

    private SplineBasis splineBasis = null;
}
