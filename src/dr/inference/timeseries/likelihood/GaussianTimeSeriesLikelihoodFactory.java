package dr.inference.timeseries.likelihood;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.timeseries.core.TimeSeriesModel;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.engine.DisabledGradientEngine;
import dr.inference.timeseries.engine.GradientEngine;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.engine.kalman.MomentAnalyticalKalmanGradientEngine;
import dr.inference.timeseries.engine.kalman.formula.CanonicalSelectionMatrixGradientFormula;
import dr.inference.timeseries.engine.kalman.CanonicalAnalyticalKalmanGradientEngine;
import dr.inference.timeseries.engine.kalman.formula.CanonicalDiffusionMatrixGradientFormula;
import dr.inference.timeseries.engine.kalman.CanonicalKalmanSmootherEngine;
import dr.inference.timeseries.engine.kalman.formula.CanonicalStationaryMeanGradientFormula;
import dr.inference.timeseries.engine.kalman.formula.CanonicalBlockDiagonalGradientCache;
import dr.inference.timeseries.engine.kalman.formula.MomentDiffusionMatrixGradientFormula;
import dr.inference.timeseries.engine.kalman.MomentSmootherResults;
import dr.inference.timeseries.engine.kalman.GaussianForwardComputationMode;
import dr.inference.timeseries.engine.kalman.GaussianLikelihoodEngineFactory;
import dr.inference.timeseries.engine.kalman.MomentKalmanSmootherEngine;
import dr.inference.timeseries.engine.kalman.formula.MomentSelectionMatrixGradientFormula;
import dr.inference.timeseries.engine.kalman.formula.MomentStationaryMeanGradientFormula;
import dr.inference.timeseries.model.gaussian.EulerOUProcessModel;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.model.gaussian.OUTimeSeriesProcessAdapter;
import dr.evomodel.continuous.ou.DiffusionMatrixParameterization;
import dr.evomodel.continuous.ou.DiffusionMatrixParameterizationFactory;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

/**
 * Top-level builder for Gaussian time-series likelihoods.
 *
 * <p>This centralizes the current policy:
 * <ul>
 *   <li>canonical forward, smoothing, and analytical gradients are the production path</li>
 *   <li>moment-form engines remain available for explicit debugging and parity checks</li>
 * </ul>
 */
public final class GaussianTimeSeriesLikelihoodFactory {

    private GaussianTimeSeriesLikelihoodFactory() {
        // no instances
    }

    public static TimeSeriesLikelihood create(final String name,
                                              final TimeSeriesModel model,
                                              final GaussianForwardComputationMode forwardMode,
                                              final GaussianGradientComputationMode gradientMode) {
        return create(name, model, forwardMode, defaultSmootherModeFor(gradientMode), gradientMode);
    }

    public static TimeSeriesLikelihood create(final String name,
                                              final TimeSeriesModel model) {
        return create(name, model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
    }

    public static TimeSeriesLikelihood create(final String name,
                                              final TimeSeriesModel model,
                                              final GaussianForwardComputationMode forwardMode,
                                              final GaussianSmootherComputationMode smootherMode,
                                              final GaussianGradientComputationMode gradientMode) {
        if (name == null) {
            throw new IllegalArgumentException("name must not be null");
        }
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        if (!(model.getObservationModel() instanceof LinearGaussianObservationModel)) {
            throw new IllegalArgumentException(
                    "Observation model must be GaussianObservationModel: " +
                            model.getObservationModel().getClass().getName());
        }

        final LatentProcessModel latentProcess = model.getLatentProcessModel();
        final RepresentableProcess process = representationFor(latentProcess);
        prepareRepeatedDeltaCache(process, model);
        final LinearGaussianObservationModel observationModel = (LinearGaussianObservationModel) model.getObservationModel();

        final LikelihoodEngine likelihoodEngine = GaussianLikelihoodEngineFactory.createForwardEngine(
                process,
                observationModel,
                model.getTimeGrid(),
                forwardMode);
        final GradientEngine gradientEngine = createGradientEngine(
                latentProcess, process, observationModel, model, smootherMode, gradientMode);
        return new TimeSeriesLikelihood(name, model, likelihoodEngine, gradientEngine);
    }

    private static RepresentableProcess representationFor(final LatentProcessModel latentProcess) {
        if (latentProcess instanceof RepresentableProcess) {
            return (RepresentableProcess) latentProcess;
        }
        throw new IllegalArgumentException(
                "Latent process must be representable for time-series inference: "
                        + latentProcess.getClass().getName());
    }

    private static GaussianSmootherComputationMode defaultSmootherModeFor(
            final GaussianGradientComputationMode gradientMode) {
        return gradientMode == GaussianGradientComputationMode.CANONICAL_ANALYTICAL
                ? GaussianSmootherComputationMode.CANONICAL
                : GaussianSmootherComputationMode.MOMENT;
    }

    private static void prepareRepeatedDeltaCache(final RepresentableProcess process,
                                                  final TimeSeriesModel model) {
        final GaussianTransitionRepresentation transitionRepresentation =
                process.getRepresentation(GaussianTransitionRepresentation.class);
        if (transitionRepresentation instanceof CachedGaussianTransitionRepresentation) {
            ((CachedGaussianTransitionRepresentation) transitionRepresentation)
                    .prepareTimeGrid(model.getTimeGrid());
        }
    }

    private static GradientEngine createGradientEngine(final LatentProcessModel latentProcess,
                                                       final RepresentableProcess process,
                                                       final LinearGaussianObservationModel observationModel,
                                                       final TimeSeriesModel model,
                                                       final GaussianSmootherComputationMode smootherMode,
                                                       final GaussianGradientComputationMode gradientMode) {
        if (gradientMode == null) {
            throw new IllegalArgumentException("gradientMode must not be null");
        }
        switch (gradientMode) {
            case DISABLED:
                return new DisabledGradientEngine();
            case MOMENT_ANALYTICAL:
                return createAnalyticalGradientEngine(latentProcess, process, observationModel, model,
                        GaussianSmootherComputationMode.MOMENT);
            case CANONICAL_ANALYTICAL:
                return createAnalyticalGradientEngine(latentProcess, process, observationModel, model, smootherMode);
            default:
                throw new IllegalArgumentException("Unsupported gradient mode: " + gradientMode);
        }
    }

    private static GradientEngine createAnalyticalGradientEngine(final LatentProcessModel latentProcess,
                                                                 final RepresentableProcess process,
                                                                 final LinearGaussianObservationModel observationModel,
                                                                 final TimeSeriesModel model,
                                                                 final GaussianSmootherComputationMode smootherMode) {
        if (latentProcess instanceof OUTimeSeriesProcessAdapter) {
            final OUProcessModel ouProcess = ((OUTimeSeriesProcessAdapter) latentProcess).getProcessModel();
            return buildAnalyticalEngineForOu(
                    ouProcess,
                    createSmoother(process, observationModel, model, smootherMode),
                    ouProcess.getDriftMatrix(),
                    ouProcess.getDiffusionMatrix(),
                    ouProcess.getStationaryMeanParameter(),
                    ouProcess.getInitialCovarianceParameter(),
                    ouProcess.getStateDimension());
        }
        if (latentProcess instanceof EulerOUProcessModel) {
            final EulerOUProcessModel eulerProcess = (EulerOUProcessModel) latentProcess;
            return buildAnalyticalEngineForOu(
                    null,
                    createSmoother(process, observationModel, model, smootherMode),
                    eulerProcess.getDriftMatrix(),
                    eulerProcess.getDiffusionMatrix(),
                    eulerProcess.getStationaryMeanParameter(),
                    eulerProcess.getInitialCovarianceParameter(),
                    eulerProcess.getStateDimension());
        }

        throw new IllegalArgumentException(
                "Moment-form analytical gradients are not configured for process type: "
                        + process.getClass().getName());
    }

    private static MomentSmootherResults createSmoother(final RepresentableProcess process,
                                                          final LinearGaussianObservationModel observationModel,
                                                          final TimeSeriesModel model,
                                                          final GaussianSmootherComputationMode smootherMode) {
        final GaussianTransitionRepresentation transitionRepresentation =
                process.getRepresentation(GaussianTransitionRepresentation.class);
        switch (smootherMode) {
            case MOMENT:
                return new MomentKalmanSmootherEngine(transitionRepresentation, observationModel, model.getTimeGrid());
            case CANONICAL:
                if (!process.supportsRepresentation(CanonicalGaussianBranchTransitionKernel.class)) {
                    throw new IllegalArgumentException(
                            "Process does not support canonical Gaussian branch transitions: "
                                    + process.getClass().getName());
                }
                return new CanonicalKalmanSmootherEngine(
                        process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                        transitionRepresentation,
                        observationModel,
                        model.getTimeGrid());
            default:
                throw new IllegalArgumentException("Unsupported smoother mode: " + smootherMode);
        }
    }

    private static GradientEngine buildAnalyticalEngineForOu(final OUProcessModel processModel,
                                                             final MomentSmootherResults smoother,
                                                             final MatrixParameterInterface driftMatrix,
                                                             final dr.inference.model.MatrixParameterInterface diffusionMatrix,
                                                             final dr.inference.model.Parameter stationaryMean,
                                                             final MatrixParameter initialCovariance,
                                                             final int stateDimension) {
        final DiffusionMatrixParameterization diffusionParameterization =
                DiffusionMatrixParameterizationFactory.create(diffusionMatrix);
        final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache =
                createBlockDiagonalGradientCache(
                        processModel,
                        driftMatrix,
                        diffusionParameterization,
                        stationaryMean,
                        initialCovariance,
                        stateDimension);
        final GradientFormulaBundle formulas = new GradientFormulaBundle(
                new MomentSelectionMatrixGradientFormula(driftMatrix, stateDimension),
                processModel != null
                        ? new CanonicalSelectionMatrixGradientFormula(
                        processModel, driftMatrix, stateDimension, blockDiagonalGradientCache)
                        : new CanonicalSelectionMatrixGradientFormula(driftMatrix, stateDimension),
                new MomentStationaryMeanGradientFormula(
                        processModel,
                        stationaryMean,
                        initialCovariance,
                        stateDimension),
                new CanonicalStationaryMeanGradientFormula(
                        processModel,
                        stationaryMean,
                        initialCovariance,
                        stateDimension,
                        blockDiagonalGradientCache),
                new MomentDiffusionMatrixGradientFormula(
                        diffusionParameterization,
                        stateDimension),
                new CanonicalDiffusionMatrixGradientFormula(
                        processModel,
                        diffusionParameterization,
                        stateDimension,
                        blockDiagonalGradientCache));
        if (smoother instanceof CanonicalKalmanSmootherEngine) {
            return new CanonicalAnalyticalKalmanGradientEngine(
                    (CanonicalKalmanSmootherEngine) smoother,
                    formulas.canonicalSelection,
                    formulas.canonicalMean,
                    formulas.canonicalDiffusion);
        }
        return new MomentAnalyticalKalmanGradientEngine(
                smoother,
                formulas.selection,
                formulas.mean,
                formulas.diffusion);
    }

    private static final class GradientFormulaBundle {
        final MomentSelectionMatrixGradientFormula selection;
        final CanonicalSelectionMatrixGradientFormula canonicalSelection;
        final MomentStationaryMeanGradientFormula mean;
        final CanonicalStationaryMeanGradientFormula canonicalMean;
        final MomentDiffusionMatrixGradientFormula diffusion;
        final CanonicalDiffusionMatrixGradientFormula canonicalDiffusion;

        GradientFormulaBundle(final MomentSelectionMatrixGradientFormula selection,
                              final CanonicalSelectionMatrixGradientFormula canonicalSelection,
                              final MomentStationaryMeanGradientFormula mean,
                              final CanonicalStationaryMeanGradientFormula canonicalMean,
                              final MomentDiffusionMatrixGradientFormula diffusion,
                              final CanonicalDiffusionMatrixGradientFormula canonicalDiffusion) {
            this.selection = selection;
            this.canonicalSelection = canonicalSelection;
            this.mean = mean;
            this.canonicalMean = canonicalMean;
            this.diffusion = diffusion;
            this.canonicalDiffusion = canonicalDiffusion;
        }
    }

    private static CanonicalBlockDiagonalGradientCache createBlockDiagonalGradientCache(
            final OUProcessModel processModel,
            final MatrixParameterInterface driftMatrix,
            final DiffusionMatrixParameterization diffusionParameterization,
            final dr.inference.model.Parameter stationaryMean,
            final MatrixParameter initialCovariance,
            final int stateDimension) {
        if (processModel == null
                || !(driftMatrix instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)
                || !CanonicalBlockDiagonalGradientCache.isAvailable(processModel, driftMatrix)) {
            return null;
        }
        return new CanonicalBlockDiagonalGradientCache(
                processModel,
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) driftMatrix,
                diffusionParameterization,
                stationaryMean,
                initialCovariance,
                stateDimension);
    }
}
