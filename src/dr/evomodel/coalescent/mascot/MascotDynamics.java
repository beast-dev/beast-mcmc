/*
 * MascotDynamics.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.model.Parameter;

/**
 * Lightweight parameter layout helper for the BEAST-X MASCOT implementation.
 */
public final class MascotDynamics {

    private final int stateCount;
    private final Parameter theta;
    private final Parameter epochTimes;

    public MascotDynamics(int stateCount, Parameter theta, Parameter epochTimes) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if (theta == null) {
            throw new IllegalArgumentException("theta parameter is required");
        }
        this.stateCount = stateCount;
        this.theta = theta;
        this.epochTimes = epochTimes;
        checkThetaDimension();
    }

    public int getStateCount() {
        return stateCount;
    }

    public Parameter getTheta() {
        return theta;
    }

    public Parameter getEpochTimes() {
        return epochTimes;
    }

    public int getEpochCount() {
        return epochTimes == null ? 1 : epochTimes.getDimension() + 1;
    }

    public int getParametersPerEpoch() {
        return stateCount * (stateCount - 1) + stateCount;
    }

    public int getParameterCount() {
        return getEpochCount() * getParametersPerEpoch();
    }

    public double[] getThetaValues() {
        checkThetaDimension();
        return theta.getParameterValues();
    }

    public double[] getBoundaries() {
        int epochCount = getEpochCount();
        double[] boundaries = new double[epochCount + 1];
        boundaries[0] = 0.0;
        if (epochTimes != null) {
            for (int i = 0; i < epochTimes.getDimension(); i++) {
                boundaries[i + 1] = epochTimes.getParameterValue(i);
            }
        }
        boundaries[boundaries.length - 1] = Double.POSITIVE_INFINITY;
        for (int i = 1; i < boundaries.length; i++) {
            if (!(boundaries[i] > boundaries[i - 1])) {
                throw new IllegalArgumentException("epoch times must be strictly increasing in backward time");
            }
        }
        return boundaries;
    }

    private void checkThetaDimension() {
        int expected = getParameterCount();
        if (theta.getDimension() != expected) {
            throw new IllegalArgumentException("theta dimension " + theta.getDimension() +
                    " does not match expected dimension " + expected);
        }
    }
}
