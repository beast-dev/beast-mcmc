/*
 * MascotDynamics.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.EpochBoundaries;
import dr.evomodel.coalescent.basta.AbstractPopulationSizeModel;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.IntervalSpecificPopulationSizeModel;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.ComplexSubstitutionModelGradientSupport;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;

import java.util.Arrays;

/**
 * Maps parser-facing migration/population parameters to GenericMascotLikelihoodDelegate's flat
 * epoch-major [migration rates, log population sizes] layout and gradients.
 */
public final class MascotDynamics {

    private final int stateCount;
    private final Parameter migrationRates;
    private final SubstitutionModel[] migrationModels;
    private final AbstractPopulationSizeModel populationSizeModel;
    private final Parameter populationSizeParameter;
    private final Parameter epochTimes;
    private final int migrationRatesPerEpoch;
    private final int parametersPerEpoch;
    private double[] migrationModelMatrix;
    private Parameter migrationModelRatesParameter;

    public MascotDynamics(int stateCount, Parameter migrationRates, AbstractPopulationSizeModel populationSizeModel,
                          Parameter epochTimes) {
        this(stateCount, migrationRates, null, populationSizeModel, epochTimes);
    }

    public MascotDynamics(int stateCount, SubstitutionModel migrationModel, AbstractPopulationSizeModel populationSizeModel,
                          Parameter epochTimes) {
        this(stateCount, null, new SubstitutionModel[]{migrationModel}, populationSizeModel, epochTimes);
    }

    public MascotDynamics(int stateCount, SubstitutionModel[] migrationModels, AbstractPopulationSizeModel populationSizeModel,
                          Parameter epochTimes) {
        this(stateCount, null, migrationModels, populationSizeModel, epochTimes);
    }

    private MascotDynamics(int stateCount, Parameter migrationRates, SubstitutionModel[] migrationModels,
                           AbstractPopulationSizeModel populationSizeModel, Parameter epochTimes) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if ((migrationRates == null) == (migrationModels == null)) {
            throw new IllegalArgumentException("exactly one migration-rate source is required: Parameter or SubstitutionModel(s)");
        }
        if (populationSizeModel == null) {
            throw new IllegalArgumentException("populationSizeModel is required");
        }
        this.stateCount = stateCount;
        this.migrationRates = migrationRates;
        this.migrationModels = migrationModels;
        this.populationSizeModel = populationSizeModel;
        this.populationSizeParameter = extractPopulationSizeParameter(populationSizeModel);
        this.epochTimes = epochTimes;
        this.migrationRatesPerEpoch = stateCount * (stateCount - 1);
        this.parametersPerEpoch = migrationRatesPerEpoch + stateCount;
        checkDimensions();
    }

    private static Parameter extractPopulationSizeParameter(AbstractPopulationSizeModel model) {
        switch (model.getModelType()) {
            case CONSTANT:
                return ((ConstantPopulationSizeModel) model).getPopulationSizeParameter();
            case PIECEWISE_CONSTANT:
                return ((IntervalSpecificPopulationSizeModel) model).getPopulationSizeParameter();
            default:
                throw new IllegalArgumentException("MascotDynamics only supports CONSTANT or " +
                        "PIECEWISE_CONSTANT population-size models, got " + model.getModelType());
        }
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
        return populationSizeParameter;
    }

    public AbstractPopulationSizeModel getPopulationSizeModel() {
        return populationSizeModel;
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

    public void writeThetaValues(double[] destination) {
        checkDimensions();
        int epochCount = getEpochCount();
        if (destination.length != epochCount * parametersPerEpoch) {
            throw new IllegalArgumentException("destination dimension " + destination.length +
                    " does not match expected dimension " + (epochCount * parametersPerEpoch));
        }
        for (int epoch = 0; epoch < epochCount; epoch++) {
            int base = epoch * parametersPerEpoch;
            if (migrationRates != null) {
                int migrationBase = epoch * migrationRatesPerEpoch;
                for (int j = 0; j < migrationRatesPerEpoch; j++) {
                    double rate = migrationRates.getParameterValue(migrationBase + j);
                    if (!(rate >= 0.0) || !Double.isFinite(rate)) {
                        throw new MascotLikelihoodDelegate.NumericalException("invalid migration rate at index " +
                                (migrationBase + j) + ": " + rate);
                    }
                    destination[base + j] = rate;
                }
            } else {
                copyMigrationModelRates(destination, base, epoch);
            }
            boolean shared = populationSizeParameter.getDimension() == stateCount;
            int popSizeBase = shared ? 0 : epoch * stateCount;
            for (int k = 0; k < stateCount; k++) {
                double naturalSize = populationSizeParameter.getParameterValue(popSizeBase + k);
                if (!(naturalSize > 0.0) || !Double.isFinite(naturalSize)) {
                    throw new MascotLikelihoodDelegate.NumericalException("invalid population size at epoch " + epoch +
                            ", state " + k + ": " + naturalSize);
                }
                destination[base + migrationRatesPerEpoch + k] = Math.log(naturalSize);
            }
        }
    }

    public double[] getThetaValues() {
        double[] theta = new double[getParameterCount()];
        writeThetaValues(theta);
        return theta;
    }

    public void writeMigrationGradient(double[] combinedGradient, double[] destination) {
        if (migrationModels != null) {
            writeMigrationModelGradient(combinedGradient, destination);
            return;
        }
        int epochCount = getEpochCount();
        if (destination.length != epochCount * migrationRatesPerEpoch) {
            throw new IllegalArgumentException("destination dimension " + destination.length +
                    " does not match expected dimension " + (epochCount * migrationRatesPerEpoch));
        }
        for (int epoch = 0; epoch < epochCount; epoch++) {
            System.arraycopy(combinedGradient, epoch * parametersPerEpoch,
                    destination, epoch * migrationRatesPerEpoch, migrationRatesPerEpoch);
        }
    }

    public double[] extractMigrationGradient(double[] combinedGradient) {
        double[] result = new double[migrationModels != null
                ? migrationModelGradientDimension() : getEpochCount() * migrationRatesPerEpoch];
        writeMigrationGradient(combinedGradient, result);
        return result;
    }

    private int migrationModelGradientDimension() {
        String compatibilityError = getMigrationGradientCompatibilityError();
        if (compatibilityError != null) {
            throw new IllegalStateException(compatibilityError);
        }
        return getMigrationModelRatesParameter().getDimension();
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

    public double[] extractPopSizeGradient(double[] combinedGradient) {
        double[] result = new double[populationSizeParameter.getDimension()];
        writePopSizeGradient(combinedGradient, result);
        return result;
    }

    public void writePopSizeGradient(double[] combinedGradient, double[] destination) {
        if (destination.length != populationSizeParameter.getDimension()) {
            throw new IllegalArgumentException("destination dimension " + destination.length +
                    " does not match expected dimension " + populationSizeParameter.getDimension());
        }
        Arrays.fill(destination, 0.0);
        int epochCount = getEpochCount();
        boolean shared = populationSizeParameter.getDimension() == stateCount;
        for (int epoch = 0; epoch < epochCount; epoch++) {
            int base = epoch * parametersPerEpoch + migrationRatesPerEpoch;
            int resultBase = shared ? 0 : epoch * stateCount;
            for (int k = 0; k < stateCount; k++) {
                double naturalSize = populationSizeParameter.getParameterValue(resultBase + k);
                destination[resultBase + k] += combinedGradient[base + k] / naturalSize;
            }
        }
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
        int expectedShared = stateCount;
        int expectedPerEpoch = epochCount * stateCount;
        int actual = populationSizeParameter.getDimension();
        if (actual != expectedShared && actual != expectedPerEpoch) {
            throw new IllegalArgumentException("population size parameter dimension " + actual +
                    " does not match expected dimension " + expectedShared + " (one shared block, via " +
                    "ConstantPopulationSizeModel) or " + expectedPerEpoch + " (one block per epoch, via " +
                    "PiecewiseConstantPopulationSizeModel)");
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
                if (!(rate >= 0.0) || !Double.isFinite(rate)) {
                    throw new MascotLikelihoodDelegate.NumericalException("invalid migrationModel rate for epoch " + epoch +
                            ", source " + source + ", sink " + sink + ": " + rate);
                }
                theta[thetaOffset + index] = rate;
                index++;
            }
        }
    }

    private void writeMigrationModelGradient(double[] combinedGradient, double[] destination) {
        String compatibilityError = getMigrationGradientCompatibilityError();
        if (compatibilityError != null) {
            throw new IllegalStateException(compatibilityError);
        }
        Parameter ratesParameter = getMigrationModelRatesParameter();
        if (destination.length != ratesParameter.getDimension()) {
            throw new IllegalArgumentException("destination dimension " + destination.length +
                    " does not match expected dimension " + ratesParameter.getDimension());
        }
        int resultOffset = 0;
        for (int epoch = 0; epoch < migrationModels.length; epoch++) {
            SubstitutionModel model = migrationModels[epoch];
            Parameter epochRatesParameter = ComplexSubstitutionModelGradientSupport.getRatesParameter(model);
            int epochGradientOffset = epoch * parametersPerEpoch;
            int parameterIndex = 0;
            for (int source = 0; source < stateCount; source++) {
                for (int sink = source + 1; sink < stateCount; sink++) {
                    destination[resultOffset + parameterIndex] = gradientWrtModelRate(
                            combinedGradient, epochGradientOffset, model, epochRatesParameter, parameterIndex, source, sink);
                    parameterIndex++;
                }
            }
            for (int sink = 0; sink < stateCount; sink++) {
                for (int source = sink + 1; source < stateCount; source++) {
                    destination[resultOffset + parameterIndex] = gradientWrtModelRate(
                            combinedGradient, epochGradientOffset, model, epochRatesParameter, parameterIndex, source, sink);
                    parameterIndex++;
                }
            }
            resultOffset += epochRatesParameter.getDimension();
        }
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
