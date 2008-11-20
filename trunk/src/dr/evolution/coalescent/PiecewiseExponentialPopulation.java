/*
 * PiecewiseLinearPopulation.java
 *
 * Copyright (C) 2002-2008 Alexei Drummond and Andrew Rambaut
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

import java.util.Collections;
import java.util.ArrayList;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewiseExponentialPopulation extends DemographicFunction.Abstract {

    /**
     * Construct demographic model with default settings
     */
    public PiecewiseExponentialPopulation(double[] intervals, double N0, double[] lambdas, Type units) {
        super(units);
        if (lambdas == null || intervals == null) {
            throw new IllegalArgumentException();
        }
        if (lambdas.length != intervals.length + 1) {
            throw new IllegalArgumentException();
        }

        this.N0 = N0;
        this.intervals = intervals;
        this.lambdas = lambdas;

    }

    // **************************************************************
    // Implementation of abstract methods
    // **************************************************************

    public double getDemographic(double t) {
        int epoch = 0;
        double t1 = t;
        double N1 = N0;
        double dt = getEpochDuration(epoch);
        while (t1 > dt) {
            N1 = N1 * Math.exp(-lambdas[epoch] * dt);
            t1 -= dt;
            epoch ++;
            dt = getEpochDuration(epoch);
        }
        return N1 * Math.exp(-lambdas[epoch] * t1);
    }

    /**
     * Gets the integral of intensity function from time 0 to time t.
     */
    public double getIntensity(double t) {
        throw new RuntimeException("Not implemented!");
    }

    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public double getIntegral(double start, double finish) {
        //final double v1 = getIntensity(finish) - getIntensity(start);
        // Until the above getIntensity is implemented, numerically integrate
        final double numerical = getNumericalIntegral(start, finish);
        return numerical;
    }

    public double getUpperBound(int i) {
        return 1e9;
    }

    public double getLowerBound(int i) {
        return Double.MIN_VALUE;
    }

    public int getNumArguments() {
        return lambdas.length + 1;
    }

    public String getArgumentName(int i) {
        if (i == 0) {
            return "N0";
        }
        return "theta" + i;
    }

    public double getArgument(int i) {
        if (i == 0) {
            return N0;
        }
        return lambdas[i - 1];
    }

    public void setArgument(int i, double value) {
        if (i == 0) {
            N0 = value;
        } else {
            lambdas[i - 1] = value;
        }
    }

    public DemographicFunction getCopy() {
        PiecewiseExponentialPopulation df = new PiecewiseExponentialPopulation(new double[intervals.length], N0, new double[lambdas.length], getUnits());
        System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
        System.arraycopy(lambdas, 0, df.lambdas, 0, lambdas.length);

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
        if (epoch < 0 || epoch >= lambdas.length) {
            throw new IllegalArgumentException();
        }
        return lambdas[epoch];
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append(N0);
        for (int i = 0; i < lambdas.length; i++) {
            buffer.append("\t").append(lambdas[i]);
        }
        return buffer.toString();
    }

    double[] intervals;
    double[] lambdas;
    double N0;
}