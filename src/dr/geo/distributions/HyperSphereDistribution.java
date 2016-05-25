/*
 * HyperSphereDistribution.java
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

import dr.math.distributions.MultivariateDistribution;
import dr.geo.math.SphericalPolarCoordinates;
import dr.geo.math.Space;

/**
 * @author Marc Suchard
 *
 * Abstract base class for probability distributions on a hyper-sphere S^{d-1} in R^{d}.
 * Usual cases are: circle (d=2) and sphere (d=3)
 */

public abstract class HyperSphereDistribution implements MultivariateDistribution {

    protected static double LOG_2_PI = Math.log(2 * Math.PI);
    
    public HyperSphereDistribution(int dim, Space space, double[] mean, double kappa) {
        this.dim = dim;
        this.mean = mean;
        this.kappa = kappa;
        this.space = space;
        checkImplementation();
    }

    public abstract double logPdf(double[] x);

    public double[][] getScaleMatrix() {
        return new double[0][];
    }

    public double[] getMean() {
        return new double[0];
    }

    public String getType() {
        return "Hyperspherical";
    }

    protected int getAllowableDim() {
        return dim;
    }

    protected void checkImplementation() {
        final int allowableDim = getAllowableDim();
        if (allowableDim != dim) {
            throw new RuntimeException(getType() + " distribution is not implemented in R^"+dim);
        }
    }


    public static double innerProduct(double[] x, double[] y, Space space) {      
        switch (space) {
            case LAT_LONG  : return latLongToCartesianInnerProduct(x, y, 1.0);
            case CARTESIAN : return cartestianInnerProduct(x, y);
            case RADIANS   : return radiansInnerProduct(x, y);
            default : throw new RuntimeException("Inner product not yet implemented");
        }
    }

    private static double radiansInnerProduct(double[] x, double[] y) {

        // x[] and y[] should be in the form (radians)
        if (x.length != 1 && y.length != 1) {
            throw new RuntimeException("Wrong dimensions");
        }

        return Math.cos(x[0] - y[0]);
    }

    private static double cartestianInnerProduct(double[] x, double[] y) {

        // x[] and y[] must be unit vectors
        if (!isUnitVector(x) || !isUnitVector(y)) {
            throw new RuntimeException("Inner produce must be on unit vectors");
        }
        
        final int len = x.length;
        double innerProduct = 0;
        for(int i = 0; i < len; i++) {
            innerProduct += x[i] * y[i];
        }
        return innerProduct;
    }

    private static double tolerance = 1E-10;

    public static boolean isUnitVector(double[] test) {
        double norm = 0;
        for(double d : test) {
            norm += d * d;
        }
        return (Math.abs(norm - 1) < tolerance);
    }

    private static double latLongToCartesianInnerProduct(double[] x, double[] y, double radius) {

        // x[] and y[] should be in the form (lat, long)
        if (x.length != 2 || y.length != 2) {
            throw new RuntimeException("Wrong dimensions");
        }

        final SphericalPolarCoordinates coordX = new SphericalPolarCoordinates(x[0], x[1], radius);
        final SphericalPolarCoordinates coordY = new SphericalPolarCoordinates(y[0], y[1], radius);

        return coordX.getCartesianCoordinates().dot(coordY.getCartesianCoordinates());
    }

    protected int dim;
    protected Space space;
    protected double[] mean;
    protected double kappa;
}
