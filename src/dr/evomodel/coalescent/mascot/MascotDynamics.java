/*
 * MascotDynamics.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.EpochBoundaries;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.ComplexSubstitutionModelGradientSupport;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;

/**
 * Lightweight parameter layout helper for the BEAST-X MASCOT implementation.
 * <p/>
 * Migration rates and population sizes are kept as separate parser-facing
 * inputs throughout (matching BASTA's {@code <popSizes>} convention).
 * Migration may be declared either as native positive rates in a {@link Parameter}
 * or through a BEAST substitution model whose realized infinitesimal matrix is
 * copied into those same natural-rate coordinates at evaluation time. {@link MascotCore}
 * still wants one flat,
 * epoch-major [migration rates, log population sizes] array per evaluation (that
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
    private final SubstitutionModel migrationModel;
    private final Parameter popSizes;
    private final Parameter epochTimes;
    private final int migrationRatesPerEpoch;
    private final int parametersPerEpoch;
    private double[] migrationModelMatrix;

    public MascotDynamics(int stateCount, Parameter migrationRates, Parameter popSizes, Parameter epochTimes) {
        this(stateCount, migrationRates, null, popSizes, epochTimes);
    }

    public MascotDynamics(int stateCount, SubstitutionModel migrationModel, Parameter popSizes, Parameter epochTimes) {
        this(stateCount, null, migrationModel, popSizes, epochTimes);
    }

    private MascotDynamics(int stateCount, Parameter migrationRates, SubstitutionModel migrationModel,
                           Parameter popSizes, Parameter epochTimes) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if ((migrationRates == null) == (migrationModel == null)) {
            throw new IllegalArgumentException("exactly one migration-rate source is required: Parameter or SubstitutionModel");
        }
        if (popSizes == null) {
            throw new IllegalArgumentException("popSizes parameter is required");
        }
        this.stateCount = stateCount;
        this.migrationRates = migrationRates;
        this.migrationModel = migrationModel;
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
        return migrationRates != null ? migrationRates : getMigrationModelRatesParameter();
    }

    public SubstitutionModel getMigrationModel() {
        return migrationModel;
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
            if (migrationRates != null) {
                int migrationBase = epoch * migrationRatesPerEpoch;
                for (int j = 0; j < migrationRatesPerEpoch; j++) {
                    double rate = migrationRates.getParameterValue(migrationBase + j);
                    if (!(rate > 0.0) || !Double.isFinite(rate)) {
                        throw new MascotCore.NumericalException("invalid migration rate at index " +
                                (migrationBase + j) + ": " + rate);
                    }
                    theta[base + j] = rate;
                }
            } else {
                copyMigrationModelRates(theta, base);
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
        if (migrationModel != null) {
            return extractMigrationModelGradient(combinedGradient);
        }
        int epochCount = getEpochCount();
        double[] result = new double[epochCount * migrationRatesPerEpoch];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            System.arraycopy(combinedGradient, epoch * parametersPerEpoch,
                    result, epoch * migrationRatesPerEpoch, migrationRatesPerEpoch);
        }
        return result;
    }

    public String getMigrationGradientCompatibilityError() {
        if (migrationRates != null) {
            return null;
        }
        return ComplexSubstitutionModelGradientSupport.getCompatibilityError(
                migrationModel, "mascotGradient part=\"migration\" with a migrationModel");
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
        return EpochBoundaries.withSentinels(epochTimes, "epochTimes");
    }

    private void checkDimensions() {
        int epochCount = getEpochCount();
        if (migrationRates != null) {
            int expectedMigration = epochCount * migrationRatesPerEpoch;
            if (migrationRates.getDimension() != expectedMigration) {
                throw new IllegalArgumentException("migration-rate (theta) parameter dimension " +
                        migrationRates.getDimension() + " does not match expected dimension " + expectedMigration);
            }
        } else {
            if (epochCount != 1) {
                throw new IllegalArgumentException("migrationModel currently supports one migration epoch; " +
                        "omit epochTimes or use parameter-backed migration rates for epoch-specific migration rates");
            }
            if (migrationModel.getDataType() == null ||
                    migrationModel.getDataType().getStateCount() != stateCount) {
                throw new IllegalArgumentException("migrationModel state count " +
                        (migrationModel.getDataType() == null ? -1 : migrationModel.getDataType().getStateCount()) +
                        " does not match mascotLikelihood stateCount " + stateCount);
            }
        }
        int expectedPopSizes = epochCount * stateCount;
        if (popSizes.getDimension() != expectedPopSizes) {
            throw new IllegalArgumentException("popSizes parameter dimension " +
                    popSizes.getDimension() + " does not match expected dimension " + expectedPopSizes);
        }
    }

    private void copyMigrationModelRates(double[] theta, int thetaOffset) {
        migrationModelMatrix = ensure(migrationModelMatrix, stateCount * stateCount);
        migrationModel.getInfinitesimalMatrix(migrationModelMatrix);
        int index = 0;
        for (int source = 0; source < stateCount; source++) {
            int row = source * stateCount;
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = migrationModelMatrix[row + sink];
                if (!(rate > 0.0) || !Double.isFinite(rate)) {
                    throw new MascotCore.NumericalException("invalid migrationModel rate for source " +
                            source + ", sink " + sink + ": " + rate);
                }
                theta[thetaOffset + index] = rate;
                index++;
            }
        }
    }

    private double[] extractMigrationModelGradient(double[] combinedGradient) {
        String compatibilityError = getMigrationGradientCompatibilityError();
        if (compatibilityError != null) {
            throw new IllegalStateException(compatibilityError);
        }
        Parameter ratesParameter = getMigrationModelRatesParameter();
        double[] result = new double[ratesParameter.getDimension()];
        int parameterIndex = 0;
        for (int source = 0; source < stateCount; source++) {
            for (int sink = source + 1; sink < stateCount; sink++) {
                result[parameterIndex] = gradientWrtModelRate(
                        combinedGradient, ratesParameter, parameterIndex, source, sink);
                parameterIndex++;
            }
        }
        for (int sink = 0; sink < stateCount; sink++) {
            for (int source = sink + 1; source < stateCount; source++) {
                result[parameterIndex] = gradientWrtModelRate(
                        combinedGradient, ratesParameter, parameterIndex, source, sink);
                parameterIndex++;
            }
        }
        return result;
    }

    private double gradientWrtModelRate(double[] combinedGradient, Parameter ratesParameter,
                                        int parameterIndex, int source, int sink) {
        double rawRate = ratesParameter.getParameterValue(parameterIndex);
        if (!(rawRate > 0.0) || !Double.isFinite(rawRate)) {
            throw new IllegalStateException("migrationModel rates parameter has a non-positive or non-finite " +
                    "entry at index " + parameterIndex + ": " + rawRate);
        }
        double scale = 1.0;
        if (((ComplexSubstitutionModel) migrationModel).getScaleRatesByFrequencies()) {
            scale = migrationModel.getFrequencyModel().getFrequency(sink);
        }
        return combinedGradient[rowMajorMigrationIndex(source, sink)] * scale;
    }

    private int rowMajorMigrationIndex(int source, int sink) {
        int index = source * (stateCount - 1) + sink;
        return sink < source ? index : index - 1;
    }

    private Parameter getMigrationModelRatesParameter() {
        return ComplexSubstitutionModelGradientSupport.getRatesParameter(migrationModel);
    }

    private static double[] ensure(double[] array, int length) {
        return array != null && array.length >= length ? array : new double[length];
    }
}
