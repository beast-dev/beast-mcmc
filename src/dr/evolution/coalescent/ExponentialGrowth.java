/*
 * ExponentialGrowth.java
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

package dr.evolution.coalescent;

/**
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 *
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class ExponentialGrowth extends ConstantPopulation {

    /**
     * Construct demographic model with default settings
     * @param units of time
     */
    public ExponentialGrowth(Type units) {

        super(units);
    }

    /**
     * @return growth rate.
     */
    public final double getGrowthRate() { return r; }

    /**
     * sets growth rate to r.
     * @param r
     */
    public void setGrowthRate(double r) { this.r = r; }

    /**
     * An alternative parameterization of this model. This
     * function sets growth rate for a given doubling time.
     * @param doublingTime
     */
    public void setDoublingTime(double doublingTime) {
        setGrowthRate( Math.log(2) / doublingTime );
    }

    // Implementation of abstract methods

    public double getDemographic(double t) {

        double r = getGrowthRate();
        if (r == 0) {
            return getN0();
        } else {
            return getN0() * Math.exp(-t * r);
        }
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    @Override
    public double getIntegral(double start, double finish) {
        double r = getGrowthRate();
        if (r == 0.0) {
            return (finish - start)/getN0();
        } else {
            // expm1 avoids cancellation when r is near zero
            return Math.exp(start*r) * Math.expm1((finish - start)*r)/getN0()/r;
        }
    }

    public double getIntensity(double t)
    {
        double r = getGrowthRate();
        if (r == 0.0) {
            return t/getN0();
        } else {
            return Math.expm1(t*r)/getN0()/r;
        }
    }

    @Override
    public double getIntensityGradient(double finishTime) {
        double r = getGrowthRate();
        if (r == 0.0) {
            return 1.0 / getN0();
        } else {
            return Math.exp(finishTime * r) / getN0();
        }
    }

    public double getIntensitySecondDerivative(double finishTime) {
        double r = getGrowthRate();
        if (r == 0.0) {
            return 0.0;
        } else {
            return Math.exp(finishTime * r) * r / getN0();
        }
    }

    @Override
    public double getLogDemographicGradient(double finishTime) {
        double r = getGrowthRate();
        if (r == 0) {
            return 0.0;
        } else {
            return -r;
        }
    }

    public double getLogDemographicSecondDerivative(double finishTime) {
        return 0;
    }

    @Override
    public double getIntensityArgumentGradient(int n, double t) {
        checkArgument(n, 2);
        if (n == 0) {
            return -getIntensity(t) / getN0();
        }
        // d/dr of (e^{rt} - 1) / (N0 r), by its Taylor series when rt is small to avoid cancellation
        final double r = getGrowthRate();
        final double rt = r * t;
        if (Math.abs(rt) < 1E-4) {
            return t * t * (0.5 + rt / 3.0 + rt * rt / 8.0) / getN0();
        }
        return (t * Math.exp(rt) / r - Math.expm1(rt) / (r * r)) / getN0();
    }

    @Override
    public double getLogDemographicArgumentGradient(int n, double t) {
        checkArgument(n, 2);
        return n == 0 ? 1.0 / getN0() : -t;
    }


    public double getInverseIntensity(double x) {

        double r = getGrowthRate();
        if (r == 0.0) {
            return getN0()*x;
        } else {
            return Math.log(1.0+getN0()*x*r)/r;
        }
    }

    public int getNumArguments() {
        return 2;
    }

    public String getArgumentName(int n) {
        if (n == 0) {
            return "N0";
        } else {
            return "r";
        }
    }

    public double getArgument(int n) {
        if (n == 0) {
            return getN0();
        } else {
            return getGrowthRate();
        }
    }

    public void setArgument(int n, double value) {
        if (n == 0) {
            setN0(value);
        } else {
            setGrowthRate(value);
        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public DemographicFunction getCopy() {
        ExponentialGrowth df = new ExponentialGrowth(getUnits());
        df.setN0(getN0());
        df.r = r;

        return df;
    }

    //
    // private stuff
    //

    private double r;
}
