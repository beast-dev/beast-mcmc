/*
 * PiecewiseExponentialPopulation.java
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

import java.util.Collections;
import java.util.ArrayList;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewiseExponentialPopulation extends DemographicFunction.Abstract {

    /**
     * Construct demographic model with default settings.
     * In this parametization there is a population size at present and a vector of
     * growth rates.
     */
    public PiecewiseExponentialPopulation(double[] intervals, double N0, double[] lambdas, Type units) {
        super(units);
        if (lambdas == null || intervals == null) {
            throw new IllegalArgumentException();
        }
        if (lambdas.length != intervals.length + 1) {
            throw new IllegalArgumentException();
        }

        this.thetas = new double[] { N0 };
        this.intervals = intervals;
        this.lambdas = lambdas;

    }

    /**
     * Construct demographic model with default settings
     * In this parametization there is a vector of population sizes and an ancestral
     * growth rate.
     */
    public PiecewiseExponentialPopulation(double[] intervals, double[] thetas, double lambda, Type units) {
        super(units);
        if (thetas == null || intervals == null) {
            throw new IllegalArgumentException();
        }
        if (thetas.length != intervals.length + 1) {
            throw new IllegalArgumentException();
        }

        this.thetas = thetas;
        this.intervals = intervals;
        this.lambdas = new double[] { lambda };

    }

    // **************************************************************
    // Implementation of abstract methods
    // **************************************************************

    public double getDemographic(double t) {
        if (thetas.length == 1) {
            int epoch = 0;
            double t1 = t;
            double N1 = thetas[0];
            double dt = getEpochDuration(epoch);
            while (t1 > dt) {
                N1 = N1 * Math.exp(-lambdas[epoch] * dt);
                t1 -= dt;
                epoch ++;
                dt = getEpochDuration(epoch);
            }
            return N1 * Math.exp(-lambdas[epoch] * t1);
        } else {
            int epoch = 0;
            double t1 = t;
            double N1 = thetas[0];
            double dt = getEpochDuration(epoch);
            while (t1 > dt) {
                epoch ++;
                N1 = thetas[epoch];
                t1 -= dt;
                dt = getEpochDuration(epoch);
            }
            double r;
            if (epoch + 1 >= thetas.length) {
                r = lambdas[0];
            } else {
                double N2 = thetas[epoch + 1];
                if (N1 < N2) {
                    r = -(Math.log(N2) - Math.log(N1)) / dt;
                } else if (N1 > N2) {
                    r = (Math.log(N1) - Math.log(N2)) / dt;
                } else {
                    r = 0.0;
                }
            }
            return N1 * Math.exp(-r * t1);
        }
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
        if (i < thetas.length) {
            return "theta" + i;
        }
        return "lambda" + (i - thetas.length);
    }

    public double getArgument(int i) {
        if (i < thetas.length) {
            return thetas[i];
        }
        return lambdas[i - thetas.length];
    }

    public void setArgument(int i, double value) {
        if (i < thetas.length) {
            thetas[i] = value;
        } else {
            lambdas[i - thetas.length] = value;
        }
    }

    public DemographicFunction getCopy() {
        PiecewiseExponentialPopulation df;

        if (thetas.length == 1) {
            df = new PiecewiseExponentialPopulation(new double[intervals.length], thetas[0], new double[lambdas.length], getUnits());
        } else {
            df = new PiecewiseExponentialPopulation(new double[intervals.length], new double[thetas.length], lambdas[0], getUnits());
        }
        System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
        System.arraycopy(thetas, 0, df.thetas, 0, thetas.length);
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

    public int getEpochCount() {
        return intervals.length + 1;
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
        if (epoch < 0 || epoch > intervals.length) {
            throw new IllegalArgumentException();
        }

        if (thetas.length == 1) {
            int ep = 0;
            double N1 = thetas[0];
            double dt = getEpochDuration(ep);
            while (epoch > ep) {
                N1 = N1 * Math.exp(-lambdas[ep] * dt);
                ep ++;
                dt = getEpochDuration(ep);
            }
            return N1;
        } else {
            return thetas[epoch];
        }
    }

    /**
     * @return the pop size of a given epoch.
     */
    public final double getEpochGrowthRate(int epoch) {
        if (epoch < 0 || epoch > intervals.length) {
            throw new IllegalArgumentException();
        }
        if (thetas.length == 1) {
            return lambdas[epoch];
        } else {
            if (epoch >= thetas.length - 1) {
                return lambdas[0];
            } else {
                double N1 = thetas[epoch];
                double N2 = thetas[epoch + 1];
                double dt = getEpochDuration(epoch);
                if (N1 < N2) {
                    return -Math.log(N2 - N1) / dt;
                } else if (N1 > N2) {
                    return Math.log(N1 - N2) / dt;
                } else {
                    return 0.0;
                }
            }
        }
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < thetas.length; i++) {
            buffer.append("\t").append(thetas[i]);
        }
        for (int i = 0; i < lambdas.length; i++) {
            buffer.append("\t").append(lambdas[i]);
        }
        return buffer.toString();
    }

    double[] intervals;
    double[] thetas;
    double[] lambdas;

}