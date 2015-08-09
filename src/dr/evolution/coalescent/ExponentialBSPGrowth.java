/*
 * ExponentialBSPGrowth.java
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
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 *
 * @version $Id: ExponentialGrowth.java,v 1.7 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class ExponentialBSPGrowth extends DemographicFunction.Abstract {

    /**
     * Construct demographic model with default settings
     * @param units of time
     */
    public ExponentialBSPGrowth(Type units) {
        super(units);
    }

    public void setup(double N0, double N1, double time) {
        this.N0 = N0;
        this.r = (Math.log(N0) - Math.log(N1)) / time;
    }

    public void setup(double N0, double r){
        this.N0 = N0;
        this.r = r;
    }

    public void setupN1(double N1, double r, double time) {
        this.r = r;
        this.N0 = N1*Math.exp(r*time);
    }

    // Implementation of abstract methods

    public double getDemographic(double t) {

        if (r == 0) {
            return N0;
        } else {
            return N0 * Math.exp(-t * r);
        }
    }

    public double getLogDemographic(double t) {
        if (r == 0) {
            return Math.log(N0);
        } else {
            return Math.log(N0) - (t * r);
        }
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    @Override
    public double getIntegral(double start, double finish) {
//        double integral1 = getNumericalIntegral(start, finish);

        double integral;
        if (r == 0.0) {
            integral = (finish - start) / N0;
        } else {
            integral = (Math.exp(finish*r) - Math.exp(start*r))/N0/r;
        }
        return integral;
    }

    /**
     * @return the number of arguments for this function.
     */
    public int getNumArguments() {
        return 0;
    }

    /**
     * @return the name of the n'th argument of this function.
     */
    public String getArgumentName(int n) {
        return null;
    }

    /**
     * @return the value of the n'th argument of this function.
     */
    public double getArgument(int n) {
        return 0;
    }

    /**
     * Sets the value of the nth argument of this function.
     */
    public void setArgument(int n, double value) {
    }

    /**
     * @return the lower bound of the nth argument of this function.
     */
    public double getLowerBound(int n) {
        return 0;
    }

    /**
     * Returns the upper bound of the nth argument of this function.
     */
    public double getUpperBound(int n) {
        return 0;
    }

    /**
     * Returns a copy of this function.
     */
    public DemographicFunction getCopy() {
        return null;
    }

    public double getIntensity(double t) {
        throw new RuntimeException("not implemented");
    }

    public double getInverseIntensity(double x) {
        throw new RuntimeException("not implemented");
    }


    private double r, N0;
}