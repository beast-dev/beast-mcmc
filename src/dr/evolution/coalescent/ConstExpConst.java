/*
 * ConstExpConst.java
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
 * This class models exponential growth from an initial population size which
 * then transitions back to a constant population size.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $ID$
 *
 */
public class ConstExpConst extends DemographicFunction.Abstract {
    public enum Parameterization {
        GROWTH_RATE,
        ANCESTRAL_POPULATION_SIZE
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstExpConst(Type units) {
        this(Parameterization.ANCESTRAL_POPULATION_SIZE, false, units);
    }

    public ConstExpConst(Parameterization parameterization, boolean useNumericalIntegrator, Type units) {
        super(units);
        this.parameterization = parameterization;
        this.useNumericalIntegrator = useNumericalIntegrator;
    }

    /**
     * @return initial population size.
     */
    public double getN0() { return N0; }

    /**
     * sets initial population size.
     * @param N0 new size
     */
    public void setN0(double N0) { this.N0 = N0; }

    /**
     * @return ancestral population size.
     */
    public double getN1() {
        if (parameterization == Parameterization.ANCESTRAL_POPULATION_SIZE) {
            return N1;
        }

        return N0 * Math.exp(-epochTime * growthRate);
    }

    /**
     * sets ancestral population size. Can only be set if the parameterization is ANCESTRAL_POPULATION_SIZE
     * @param N1 new size
     */
    public void setN1(double N1) {
        assert(parameterization == Parameterization.ANCESTRAL_POPULATION_SIZE);
        this.N1 = N1;
    }

    /**
     * return the first transition time to exponential growth
     * @return
     */
    public double getTime1() {
        return time1;
    }

    public void setTime1(double time1) {
        this.time1 = time1;
    }

    /**
     * Get the time period of exponential growth phase
     * @return
     */
    public double getEpochTime() {
        return epochTime;
    }

    public void setEpochTime(double epochTime) {
        this.epochTime = epochTime;
    }

    /**
     * Return the second transition time to constant size N1
     * @return
     */
    public double getTime2() {
        return time1 + epochTime;
    }

    /**
     * sets growth rate for exponential phase. Can only be set if the parameterization is GROWTH_RATE
     * @param growthRate new growthRate
     */
    public final void setGrowthRate(double growthRate) {
        assert(parameterization == Parameterization.GROWTH_RATE);
        this.growthRate = growthRate;
    }

    /**
     * @return growth rate.
     */
    public final double getGrowthRate() {
        if (parameterization == Parameterization.GROWTH_RATE) {
            return growthRate;
        }

        return (Math.log(N0) - Math.log(N1)) / epochTime;
    }


    // Implementation of abstract methods

    public double getDemographic(double t) {

        double r = getGrowthRate();

        if (t < getTime1()) {
            return getN0();
        }

        if (t >= getTime2()) {
            return getN1();
        }

        return getN0() * Math.exp(-r*(t - getTime1()));
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    public double getIntegral(double start, double finish)
    {
        if (useNumericalIntegrator) {
            return getNumericalIntegral(start, finish);
        } else {
            return getIntensity(finish) - getIntensity(start);
        }
    }

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        double time2 = getTime2();
        double oneOverN0 = 1.0 / getN0();
        double r = getGrowthRate();

        if (t < time1) {
            return (t * oneOverN0);
        }
        if (t > time1 && t < time2) {
            return (time1 * oneOverN0) + (( (Math.exp(t*r) - Math.exp(time1*r)) * oneOverN0) / r);
        }

        double oneOverN1 = 1.0 / getN1();
        // if (t >= time2) {
        return (time1 * oneOverN0) + (( (Math.exp(time2*r) - Math.exp(time1*r)) * oneOverN0) / r) + (oneOverN1 * (t-time2));
    }

    public double getInverseIntensity(double x) {

        throw new RuntimeException("Not implemented!");
    }

    public int getNumArguments() {
        return 4;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0: return "N0";
            case 1: return (parameterization == Parameterization.ANCESTRAL_POPULATION_SIZE ? "N1": "r");
            case 2: return "epochTime";
            case 3: return "time1";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0: return getN0();
            case 1: return (parameterization == Parameterization.ANCESTRAL_POPULATION_SIZE ? getN1(): getGrowthRate());
            case 2: return getEpochTime();
            case 3: return getTime1();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0: setN0(value); break;
            case 1:
                if (parameterization == Parameterization.ANCESTRAL_POPULATION_SIZE) {
                    setN1(value);
                } else {
                    setGrowthRate(value);
                }
                break;
            case 2: setEpochTime(value); break;
            case 3: setTime1(value); break;
            default: throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        if (n == 1 && parameterization == Parameterization.GROWTH_RATE) {
            return Double.NEGATIVE_INFINITY;
        }
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    //
    // private stuff
    //

    private double N0;
    private double N1;
    private double time1;
    private double epochTime;
    private double growthRate;
    private final boolean useNumericalIntegrator;
    private final Parameterization parameterization;
}