/*
 * PiecewiseConstantPopulation.java
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

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewiseConstantPopulation extends DemographicFunction.Abstract {

    ArrayList<Double> cif = null;
    ArrayList<Double> endTime = null;
    final boolean cacheCumulativeIntensities = true;

    /**
     * Construct demographic model with default settings
     *
     * @param units of time
     */
    public PiecewiseConstantPopulation(Type units) {
        super(units);
    }

    /**
     * Creates a piecewise constant model with the given break points.
     *
     * @param intervals an array of intervals, Each interval represents time
     *                  over which a different population size is assumed, the first interval starts from the
     *                  present and proceeds back into the past.
     * @param thetas    the population sizes of each interval. The last population size represents the
     *                  time from the end of the last interval out to infinity.
     */
    public PiecewiseConstantPopulation(double[] intervals, double[] thetas, Type units) {
        super(units);
        setIntervals(intervals, thetas);
    }

    public void setIntervals(final double[] intervals, final double[] thetas) {
        if (thetas == null || intervals == null) {
            throw new IllegalArgumentException();
        }
        if (thetas.length != intervals.length + 1) {
            throw new IllegalArgumentException();
        }

        this.intervals = intervals;
        this.thetas = thetas;

        if (cacheCumulativeIntensities) {
            cif = new ArrayList<Double>(intervals.length);
            endTime = new ArrayList<Double>(intervals.length);
            cif.add(getIntensity(0));
            endTime.add(intervals[0]);
            for (int i = 1; i < intervals.length; i++) {
                cif.add(getIntensity(i) + cif.get(i - 1));
                endTime.add(intervals[i] + endTime.get(i - 1));
            }
        }
    }

    // **************************************************************
    // Implementation of abstract methods
    // **************************************************************

    public double getDemographic(double t) {
        int epoch = 0;
        double t1 = t;
        while (t1 > getEpochDuration(epoch)) {
            t1 -= getEpochDuration(epoch);
            epoch += 1;
        }
        return getDemographic(epoch, t1);
    }

    /**
     * Gets the integral of intensity function from time 0 to time t.
     */
    public double getIntensity(double t) {

        if (cif != null) {

            int epoch = Collections.binarySearch(endTime, t);

            if (epoch < 0) {
                epoch = -epoch - 1;

                if (epoch > 0) {

                    return cif.get(epoch - 1) + getIntensity(epoch, t - endTime.get(epoch - 1));
                } else {
                    assert epoch == 0;
                    return getIntensity(0, t);
                }
            } else {
                return cif.get(epoch);
            }
        } else {

            double intensity = 0.0;
            int epoch = 0;
            double t1 = t;
            while (t1 > getEpochDuration(epoch)) {
                t1 -= getEpochDuration(epoch);
                intensity += getIntensity(epoch);
                epoch += 1;
            }

            // probably same bug as in EmpiricalPiecewiseConstant when t1 goes negative
            // add last fraction of intensity
            intensity += getIntensity(epoch, t1);

            return intensity;
        }
    }

    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public double getUpperBound(int i) {
        return 1e9;
    }

    public double getLowerBound(int i) {
        return Double.MIN_VALUE;
    }

    public int getNumArguments() {
        return thetas.length;
    }

    public String getArgumentName(int i) {
        return "theta" + i;
    }

    public double getArgument(int i) {
        return thetas[i];
    }

    public void setArgument(int i, double value) {
        thetas[i] = value;
    }

    public DemographicFunction getCopy() {
        PiecewiseConstantPopulation df = new PiecewiseConstantPopulation(new double[intervals.length], new double[thetas.length], getUnits());
        System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
        System.arraycopy(thetas, 0, df.thetas, 0, thetas.length);

        return df;
    }


    /**
     * @return the value of the demographic function for the given epoch and time relative to start of epoch.
     */
    protected double getDemographic(int epoch, double t) {
        return getEpochDemographic(epoch);
    }

    /**
     * @return the value of the intensity function for the given epoch.
     */
    protected double getIntensity(int epoch) {
        return getEpochDuration(epoch) / getEpochDemographic(epoch);
    }

    /**
     * @return the value of the intensity function for the given epoch and time relative to start of epoch.
     */
    protected double getIntensity(final int epoch, final double relativeTime) {
        assert 0 <= relativeTime && relativeTime <= getEpochDuration(epoch);

        return relativeTime / getEpochDemographic(epoch);
    }

    /**
     * @return the duration of the specified epoch (in whatever units this demographic model is specified in).
     */
    public double getEpochDuration(int epoch) {
        if (epoch == intervals.length) {
            return Double.POSITIVE_INFINITY;
        } else return intervals[epoch];
    }

    public void setEpochDuration(int epoch, double duration) {
        if (epoch < 0 || epoch >= intervals.length) {
            throw new IllegalArgumentException("epoch must be between 0 and " + (intervals.length - 1));
        }
        if (duration < 0.0) {
            throw new IllegalArgumentException("duration must be positive.");
        }
        intervals[epoch] = duration;
    }

    /**
     * @return the pop size of a given epoch.
     */
    public final double getEpochDemographic(int epoch) {
        if (epoch < 0 || epoch >= thetas.length) {
            throw new IllegalArgumentException();
        }
        return thetas[epoch];
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append(thetas[0]);
        for (int i = 1; i < thetas.length; i++) {
            buffer.append("\t").append(thetas[i]);
        }
        return buffer.toString();
    }

    double[] intervals;
    double[] thetas;
}
