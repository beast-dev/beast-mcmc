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
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;

/**
 * Lightweight parameter layout helper for the BEAST-X MASCOT implementation.
 * <p/>
 * Migration rates and population sizes are kept as separate parser-facing
 * inputs throughout (matching BASTA's {@code <popSizes>} convention).
 * Migration may be declared either as native positive rates in a {@link Parameter}
 * (one flat, epoch-major-concatenated array covering every epoch) or through
 * one {@link SubstitutionModel} per epoch, each supplying that epoch's
 * realized infinitesimal matrix. {@link MascotCore} still wants one flat,
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
    // One per epoch when migrationRates == null; null when migrationRates != null.
    private final SubstitutionModel[] migrationModels;
    private final Parameter popSizes;
    private final Parameter epochTimes;
    private final int migrationRatesPerEpoch;
    private final int parametersPerEpoch;
    private double[] migrationModelMatrix;
    // Lazily built, cached: a stable Parameter identity is required across
    // repeated getMigrationRates() calls (HMC operators hold onto the
    // returned reference for the whole analysis), and migrationModels.length
    // == 1 returns that model's own rates parameter directly rather than
    // wrapping a single element for no reason.
    private Parameter migrationModelRatesParameter;

    public MascotDynamics(int stateCount, Parameter migrationRates, Parameter popSizes, Parameter epochTimes) {
        this(stateCount, migrationRates, null, popSizes, epochTimes);
    }

    public MascotDynamics(int stateCount, SubstitutionModel migrationModel, Parameter popSizes, Parameter epochTimes) {
        this(stateCount, null, new SubstitutionModel[]{migrationModel}, popSizes, epochTimes);
    }

    public MascotDynamics(int stateCount, SubstitutionModel[] migrationModels, Parameter popSizes, Parameter epochTimes) {
        this(stateCount, null, migrationModels, popSizes, epochTimes);
    }

    private MascotDynamics(int stateCount, Parameter migrationRates, SubstitutionModel[] migrationModels,
                           Parameter popSizes, Parameter epochTimes) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if ((migrationRates == null) == (migrationModels == null)) {
            throw new IllegalArgumentException("exactly one migration-rate source is required: Parameter or SubstitutionModel(s)");
        }
        if (popSizes == null) {
            throw new IllegalArgumentException("popSizes parameter is required");
        }
        this.stateCount = stateCount;
        this.migrationRates = migrationRates;
        this.migrationModels = migrationModels;
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

    public SubstitutionModel[] getMigrationModels() {
        return migrationModels;
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
                copyMigrationModelRates(theta, base, epoch);
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
        if (migrationModels != null) {
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
        for (int epoch = 0; epoch < migrationModels.length; epoch++) {
            String error = ComplexSubstitutionModelGradientSupport.getCompatibilityError(
                    migrationModels[epoch], "mascotGradient part=\"migration\" with epoch " + epoch + "'s migrationModel");
            if (error != null) {
                return error;
            }
        }
        // A CompoundParameter's dimensions must each correspond to exactly one
        // underlying value slot; two epochs sharing one rates Parameter would
        // break that 1:1 mapping (and would mean an HMC move on one epoch's
        // "copy" silently also moves the other's).
        for (int i = 0; i < migrationModels.length; i++) {
            Parameter ratesI = ComplexSubstitutionModelGradientSupport.getRatesParameter(migrationModels[i]);
            for (int j = i + 1; j < migrationModels.length; j++) {
                Parameter ratesJ = ComplexSubstitutionModelGradientSupport.getRatesParameter(migrationModels[j]);
                if (ratesI == ratesJ) {
                    return "mascotGradient part=\"migration\" requires a distinct rates Parameter per epoch's " +
                            "migrationModel; epochs " + i + " and " + j + " share one";
                }
            }
        }
        return null;
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
            if (migrationModels.length != epochCount) {
                throw new IllegalArgumentException("migrationModel count (" + migrationModels.length +
                        ") does not match epoch count (" + epochCount + "); supply exactly one migrationModel " +
                        "per epoch (a single model when gridPoints is omitted)");
            }
            for (int epoch = 0; epoch < migrationModels.length; epoch++) {
                SubstitutionModel model = migrationModels[epoch];
                if (model.getDataType() == null || model.getDataType().getStateCount() != stateCount) {
                    throw new IllegalArgumentException("migrationModel for epoch " + epoch + " has state count " +
                            (model.getDataType() == null ? -1 : model.getDataType().getStateCount()) +
                            ", which does not match mascotLikelihood stateCount " + stateCount);
                }
            }
        }
        int expectedPopSizes = epochCount * stateCount;
        if (popSizes.getDimension() != expectedPopSizes) {
            throw new IllegalArgumentException("popSizes parameter dimension " +
                    popSizes.getDimension() + " does not match expected dimension " + expectedPopSizes);
        }
    }

    private void copyMigrationModelRates(double[] theta, int thetaOffset, int epoch) {
        migrationModelMatrix = ensure(migrationModelMatrix, stateCount * stateCount);
        migrationModels[epoch].getInfinitesimalMatrix(migrationModelMatrix);
        int index = 0;
        for (int source = 0; source < stateCount; source++) {
            int row = source * stateCount;
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = migrationModelMatrix[row + sink];
                if (!(rate > 0.0) || !Double.isFinite(rate)) {
                    throw new MascotCore.NumericalException("invalid migrationModel rate for epoch " + epoch +
                            ", source " + source + ", sink " + sink + ": " + rate);
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
        int resultOffset = 0;
        for (int epoch = 0; epoch < migrationModels.length; epoch++) {
            SubstitutionModel model = migrationModels[epoch];
            Parameter epochRatesParameter = ComplexSubstitutionModelGradientSupport.getRatesParameter(model);
            int epochGradientOffset = epoch * parametersPerEpoch;
            int parameterIndex = 0;
            for (int source = 0; source < stateCount; source++) {
                for (int sink = source + 1; sink < stateCount; sink++) {
                    result[resultOffset + parameterIndex] = gradientWrtModelRate(
                            combinedGradient, epochGradientOffset, model, epochRatesParameter, parameterIndex, source, sink);
                    parameterIndex++;
                }
            }
            for (int sink = 0; sink < stateCount; sink++) {
                for (int source = sink + 1; source < stateCount; source++) {
                    result[resultOffset + parameterIndex] = gradientWrtModelRate(
                            combinedGradient, epochGradientOffset, model, epochRatesParameter, parameterIndex, source, sink);
                    parameterIndex++;
                }
            }
            resultOffset += epochRatesParameter.getDimension();
        }
        return result;
    }

    private double gradientWrtModelRate(double[] combinedGradient, int epochGradientOffset, SubstitutionModel model,
                                        Parameter ratesParameter, int parameterIndex, int source, int sink) {
        double rawRate = ratesParameter.getParameterValue(parameterIndex);
        if (!(rawRate > 0.0) || !Double.isFinite(rawRate)) {
            throw new IllegalStateException("migrationModel rates parameter has a non-positive or non-finite " +
                    "entry at index " + parameterIndex + ": " + rawRate);
        }
        double scale = 1.0;
        if (((ComplexSubstitutionModel) model).getScaleRatesByFrequencies()) {
            scale = model.getFrequencyModel().getFrequency(sink);
        }
        return combinedGradient[epochGradientOffset + rowMajorMigrationIndex(source, sink)] * scale;
    }

    private int rowMajorMigrationIndex(int source, int sink) {
        int index = source * (stateCount - 1) + sink;
        return sink < source ? index : index - 1;
    }

    private Parameter getMigrationModelRatesParameter() {
        if (migrationModels.length == 1) {
            return ComplexSubstitutionModelGradientSupport.getRatesParameter(migrationModels[0]);
        }
        if (migrationModelRatesParameter == null) {
            Parameter[] ratesParameters = new Parameter[migrationModels.length];
            for (int epoch = 0; epoch < migrationModels.length; epoch++) {
                Parameter rates = ComplexSubstitutionModelGradientSupport.getRatesParameter(migrationModels[epoch]);
                if (rates == null) {
                    // Let getMigrationGradientCompatibilityError() report the
                    // actual reason at the call site instead of NPE-ing here.
                    return null;
                }
                ratesParameters[epoch] = rates;
            }
            migrationModelRatesParameter = new CompoundParameter("migrationModelRates", ratesParameters);
        }
        return migrationModelRatesParameter;
    }

    private static double[] ensure(double[] array, int length) {
        return array != null && array.length >= length ? array : new double[length];
    }
}
