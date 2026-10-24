/*
 * Expansion.java
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
 * This class models exponential growth from an initial ancestral population size.
 * (Parameters: N0=present-day population size; N1=ancestral population size; r=growth rate).
 * This model is nested with the exponential-growth population size model (N1=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class Expansion extends ExponentialGrowth {

    /**
     * Construct demographic model with default settings
     *
     * @param units of time
     */
    public Expansion(Type units) {

        super(units);
    }

    public double getN1() {
        return N1;
    }

    public void setN1(double N1) {
        this.N1 = N1;
    }

    public void setProportion(double p) {
        this.N1 = getN0() * p;
    }

    // Implementation of abstract methods

    public double getDemographic(double t) {

        double N0 = getN0();
        double N1 = getN1();
        double r = getGrowthRate();

        if (N1 > N0) throw new IllegalArgumentException("N0 must be greater than N1!");

        return N1 + ((N0 - N1) * Math.exp(-r * t));
    }

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        double N0 = getN0();
        double N1 = getN1();
        double b = (N0 - N1);
        double r = getGrowthRate();

        return Math.log(b + N1 * Math.exp(r * t)) / (r * N1);
    }

    // N(t) = N1 + (N0 - N1) e^{-rt}; the exponential-growth overrides do not apply when N1 > 0
    @Override
    public double getIntensityGradient(double finishTime) {
        return 1.0 / getDemographic(finishTime);
    }

    @Override
    public double getIntensitySecondDerivative(double finishTime) {
        final double demographic = getDemographic(finishTime);
        return getGrowthRate() * decayingPart(finishTime) / (demographic * demographic);
    }

    @Override
    public double getLogDemographicGradient(double finishTime) {
        return -getGrowthRate() * decayingPart(finishTime) / getDemographic(finishTime);
    }

    @Override
    public double getLogDemographicSecondDerivative(double finishTime) {
        final double r = getGrowthRate();
        final double demographic = getDemographic(finishTime);
        return r * r * decayingPart(finishTime) * getN1() / (demographic * demographic);
    }

    // I(t) = log D / (r N1) with D = N0 - N1 + N1 e^{rt}
    @Override
    public double getIntensityArgumentGradient(int n, double t) {
        checkArgument(n, 3);
        final double N1 = getN1();
        final double r = getGrowthRate();
        final double ert = Math.exp(r * t);
        final double D = getN0() - N1 + N1 * ert;
        switch (n) {
            case 0:
                return 1.0 / (D * r * N1);
            case 1:
                return t * ert / (r * D) - Math.log(D) / (r * r * N1);
            default:
                return Math.expm1(r * t) / (D * r * N1) - Math.log(D) / (r * N1 * N1);
        }
    }

    @Override
    public double getLogDemographicArgumentGradient(int n, double t) {
        checkArgument(n, 3);
        final double r = getGrowthRate();
        final double demographic = getDemographic(t);
        switch (n) {
            case 0:
                return Math.exp(-r * t) / demographic;
            case 1:
                return -t * decayingPart(t) / demographic;
            default:
                return -Math.expm1(-r * t) / demographic;
        }
    }

    private double decayingPart(double t) {
        return (getN0() - getN1()) * Math.exp(-getGrowthRate() * t);
    }

    public double getInverseIntensity(double x) {

        /* AER - I think this is right but until someone checks it...
          double nZero = getN0();
          double nOne = getN1();
          double r = getGrowthRate();

          if (r == 0) {
              return nZero*x;
          } else if (alpha == 0) {
              return Math.log(1.0+nZero*x*r)/r;
          } else {
              return Math.log(-(nOne/nZero) + Math.exp(nOne*x*r))/r;
          }
          */
        throw new RuntimeException("Not implemented!");
    }

    public double getIntegral(double start, double finish) {
        double v1 = getIntensity(finish) - getIntensity(start);
        //double v1 =  getNumericalIntegral(start, finish);

        return v1;
    }

    public int getNumArguments() {
        return 3;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N0";
            case 1:
                return "r";
            case 2:
                return "N1";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getN1();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0:
                setN0(value);
                break;
            case 1:
                setGrowthRate(value);
                break;
            case 2:
                setN1(value);
                break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    //
    // private stuff
    //

    private double N1 = 0.0;
}
