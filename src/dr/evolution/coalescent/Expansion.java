/*
 * Expansion.java
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

package dr.evolution.coalescent;

/**
 * This class models exponential growth from an initial ancestral population size.
 * (Parameters: N0=present-day population size; N1=ancestral population size; r=growth rate).
 * This model is nested with the exponential-growth population size model (N1=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: Expansion.java,v 1.4 2005/05/24 20:25:56 rambaut Exp $
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
