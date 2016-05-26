/*
 * PiecewiseLinearPopulation.java
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

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PiecewiseLinearPopulation extends PiecewiseConstantPopulation {

    /**
     * Construct demographic model with default settings
     *
     * @param units of time
     */
    public PiecewiseLinearPopulation(Type units) {
        super(units);
    }

    /**
     * Construct demographic model with default settings
     */
    public PiecewiseLinearPopulation(double[] intervals, double[] thetas, Type units) {

        super(intervals, thetas, units);
    }

    // **************************************************************
    // Implementation of abstract methods
    // **************************************************************

    /**
     * @return the value of the demographic function for the given epoch and time relative to start of epoch.
     */
    protected final double getDemographic(int epoch, double t) {

        // if in last epoch then the population is flat.
        if (epoch == (thetas.length - 1)) {
            return getEpochDemographic(epoch);
        }

        final double popSize1 = getEpochDemographic(epoch);
        final double popSize2 = getEpochDemographic(epoch + 1);

        final double width = getEpochDuration(epoch);

        assert 0 <= t && t <= width;

        return popSize1 + (t / width) * (popSize2 - popSize1);
    }

    public DemographicFunction getCopy() {
        PiecewiseLinearPopulation df = new PiecewiseLinearPopulation(new double[intervals.length], new double[thetas.length], getUnits());
        System.arraycopy(intervals, 0, df.intervals, 0, intervals.length);
        System.arraycopy(thetas, 0, df.thetas, 0, thetas.length);

        return df;
    }

    /**
     * @return the value of the intensity function for the given epoch.
     */
    protected final double getIntensity(int epoch) {

        final double N1 = getEpochDemographic(epoch);
        final double N2 = getEpochDemographic(epoch + 1);
        final double w = getEpochDuration(epoch);

        if (N1 != N2) {
            return w * Math.log(N2 / N1) / (N2 - N1);
        } else {
            return w / N1;
        }
    }

    /**
     * @return the value of the intensity function for the given epoch and time relative to start of epoch.
     */
    protected final double getIntensity(int epoch, double relativeTime) {
        assert relativeTime <= getEpochDuration(epoch);

        final double N1 = getEpochDemographic(epoch);
        final double N2 = getEpochDemographic(epoch + 1);
        final double w = getEpochDuration(epoch);

        if (N1 != N2) {
            return w * Math.log(N1 * w / (N2 * relativeTime + N1 * (w - relativeTime))) / (N1 - N2);
        } else {
            return relativeTime / N1;
        }
    }

    /**
     * @param targetI the intensity
     * @return the time corresponding to the given target intensity
     */
    public final double getInverseIntensity(double targetI) {

        if (cif != null) {

            int epoch = Collections.binarySearch(cif, targetI);

            if (epoch < 0) {
                epoch = -epoch - 1;

                if (epoch > 0) {

                    return endTime.get(epoch - 1) +
                            getInverseIntensity(epoch, targetI - cif.get(epoch - 1));
                } else {
                    assert epoch == 0;
                    return getInverseIntensity(0, targetI);
                }
            } else {
                return endTime.get(epoch);
            }
        } else {

            int epoch = 0;
            double cI = 0;
            double eI = getIntensity(epoch);
            double time = 0.0;
            while (cI + eI < targetI) {
                cI += eI;
                time += getEpochDuration(epoch);

                epoch += 1;
                eI = getIntensity(epoch);
            }

            time += getInverseIntensity(epoch, targetI - cI);

            return time;
        }
    }

    /**
     * @param epoch the epoch
     * @param I     intensity corresponding to relative time within this epoch
     * @return the relative time within the given epoch to accumulate the given residual intensity
     */
    private double getInverseIntensity(int epoch, double I) {

        final double N1 = getEpochDemographic(epoch);
        final double N2 = getEpochDemographic(epoch + 1);
        final double w = getEpochDuration(epoch);

        final double dn = N2 - N1;
        double time;
        if (dn != 0.0) {
            time = N1 * w * (Math.exp(I * dn / w) - 1.0) / dn;
        } else {
            time = I * N1;
        }

        return time;
    }
}
