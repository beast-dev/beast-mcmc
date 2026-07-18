/*
 * MascotDynamics.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.model.Parameter;

/**
 * Lightweight parameter layout helper for the BEAST-X MASCOT implementation.
 * <p/>
 * Migration rates and population sizes are kept as two separate parameters
 * throughout (matching BASTA's {@code <popSizes>} convention), each laid out
 * epoch-major within itself. {@link MascotCore} still wants one flat,
 * epoch-major [migration rates, population sizes] array per evaluation (that
 * layout is baked into its forward ODE and reverse-mode adjoint); this class
 * interleaves the two source parameters into that layout on the way in
 * ({@link #getThetaValues()}) and de-interleaves {@link MascotCore}'s
 * returned combined gradient back into per-parameter slices on the way out
 * ({@link #extractMigrationGradient}/{@link #extractPopSizeGradient}), so
 * migration rates and population sizes can be exposed as independent
 * {@link dr.inference.hmc.GradientWrtParameterProvider}s (see
 * {@link MascotGradient}) without MascotCore's own math changing at all.
 */
public final class MascotDynamics {

    private final int stateCount;
    private final Parameter migrationRates;
    private final Parameter popSizes;
    private final Parameter epochTimes;
    private final int migrationRatesPerEpoch;
    private final int parametersPerEpoch;

    public MascotDynamics(int stateCount, Parameter migrationRates, Parameter popSizes, Parameter epochTimes) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if (migrationRates == null) {
            throw new IllegalArgumentException("migration-rate (theta) parameter is required");
        }
        if (popSizes == null) {
            throw new IllegalArgumentException("popSizes parameter is required");
        }
        this.stateCount = stateCount;
        this.migrationRates = migrationRates;
        this.popSizes = popSizes;
        this.epochTimes = epochTimes;
        this.migrationRatesPerEpoch = stateCount * (stateCount - 1);
        this.parametersPerEpoch = migrationRatesPerEpoch + stateCount;
        checkDimensions();
    }

    public int getStateCount() {
        return stateCount;
    }

    public Parameter getMigrationRates() {
        return migrationRates;
    }

    public Parameter getPopSizes() {
        return popSizes;
    }

    public Parameter getEpochTimes() {
        return epochTimes;
    }

    public int getEpochCount() {
        return epochTimes == null ? 1 : epochTimes.getDimension() + 1;
    }

    public int getMigrationRatesPerEpoch() {
        return migrationRatesPerEpoch;
    }

    public int getParametersPerEpoch() {
        return parametersPerEpoch;
    }

    public int getParameterCount() {
        return getEpochCount() * parametersPerEpoch;
    }

    /**
     * Interleaves migrationRates/popSizes into MascotCore's expected flat,
     * epoch-major [migration rates, population sizes] layout.
     */
    public double[] getThetaValues() {
        checkDimensions();
        int epochCount = getEpochCount();
        double[] theta = new double[epochCount * parametersPerEpoch];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            int base = epoch * parametersPerEpoch;
            int migrationBase = epoch * migrationRatesPerEpoch;
            for (int j = 0; j < migrationRatesPerEpoch; j++) {
                theta[base + j] = migrationRates.getParameterValue(migrationBase + j);
            }
            int popSizeBase = epoch * stateCount;
            for (int k = 0; k < stateCount; k++) {
                theta[base + migrationRatesPerEpoch + k] = popSizes.getParameterValue(popSizeBase + k);
            }
        }
        return theta;
    }

    /**
     * De-interleaves a flat, epoch-major combined gradient (MascotCore's own
     * output layout) into just the migration-rate slice, in migrationRates'
     * own (epoch-major-concatenated) layout.
     */
    public double[] extractMigrationGradient(double[] combinedGradient) {
        int epochCount = getEpochCount();
        double[] result = new double[epochCount * migrationRatesPerEpoch];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            System.arraycopy(combinedGradient, epoch * parametersPerEpoch,
                    result, epoch * migrationRatesPerEpoch, migrationRatesPerEpoch);
        }
        return result;
    }

    /** Same as {@link #extractMigrationGradient} but for the population-size slice. */
    public double[] extractPopSizeGradient(double[] combinedGradient) {
        int epochCount = getEpochCount();
        double[] result = new double[epochCount * stateCount];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            System.arraycopy(combinedGradient, epoch * parametersPerEpoch + migrationRatesPerEpoch,
                    result, epoch * stateCount, stateCount);
        }
        return result;
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

    private void checkDimensions() {
        int epochCount = getEpochCount();
        int expectedMigration = epochCount * migrationRatesPerEpoch;
        if (migrationRates.getDimension() != expectedMigration) {
            throw new IllegalArgumentException("migration-rate (theta) parameter dimension " +
                    migrationRates.getDimension() + " does not match expected dimension " + expectedMigration);
        }
        int expectedPopSizes = epochCount * stateCount;
        if (popSizes.getDimension() != expectedPopSizes) {
            throw new IllegalArgumentException("popSizes parameter dimension " +
                    popSizes.getDimension() + " does not match expected dimension " + expectedPopSizes);
        }
    }
}
