package dr.inference.timeseries.likelihood;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.timeseries.core.TimeSeriesModel;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.engine.DisabledGradientEngine;
import dr.inference.timeseries.engine.GradientEngine;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.engine.gaussian.AnalyticalKalmanGradientEngine;
import dr.inference.timeseries.engine.gaussian.CanonicalSelectionMatrixGradientFormula;
import dr.inference.timeseries.engine.gaussian.CanonicalAnalyticalKalmanGradientEngine;
import dr.inference.timeseries.engine.gaussian.CanonicalDiffusionMatrixGradientFormula;
import dr.inference.timeseries.engine.gaussian.CanonicalKalmanSmootherEngine;
import dr.inference.timeseries.engine.gaussian.CanonicalStationaryMeanGradientFormula;
import dr.inference.timeseries.engine.gaussian.DiffusionMatrixGradientFormula;
import dr.inference.timeseries.engine.gaussian.GaussianSmootherResults;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.engine.gaussian.GaussianLikelihoodEngineFactory;
import dr.inference.timeseries.engine.gaussian.KalmanSmootherEngine;
import dr.inference.timeseries.engine.gaussian.SelectionMatrixGradientFormula;
import dr.inference.timeseries.engine.gaussian.StationaryMeanGradientFormula;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterizationFactory;
import dr.inference.timeseries.gaussian.EulerOUProcessModel;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

/**
 * Top-level builder for Gaussian time-series likelihoods.
 *
 * <p>This centralizes the current policy:
 * <ul>
 *   <li>forward likelihoods can be built in expectation or canonical form</li>
 *   <li>analytical gradients, when requested, are still computed on the established
 *       expectation/smoother backend</li>
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
        return create(name, model, forwardMode, GaussianSmootherComputationMode.EXPECTATION, gradientMode);
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
        if (!(model.getObservationModel() instanceof GaussianObservationModel)) {
            throw new IllegalArgumentException(
                    "Observation model must be GaussianObservationModel: " +
                            model.getObservationModel().getClass().getName());
        }

        final LatentProcessModel latentProcess = model.getLatentProcessModel();
        final RepresentableProcess process = representationFor(latentProcess);
        final GaussianObservationModel observationModel = (GaussianObservationModel) model.getObservationModel();

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

    private static GradientEngine createGradientEngine(final LatentProcessModel latentProcess,
                                                       final RepresentableProcess process,
                                                       final GaussianObservationModel observationModel,
                                                       final TimeSeriesModel model,
                                                       final GaussianSmootherComputationMode smootherMode,
                                                       final GaussianGradientComputationMode gradientMode) {
        if (gradientMode == null) {
            throw new IllegalArgumentException("gradientMode must not be null");
        }
        switch (gradientMode) {
            case DISABLED:
                return new DisabledGradientEngine();
            case EXPECTATION_ANALYTICAL:
                return createAnalyticalGradientEngine(latentProcess, process, observationModel, model,
                        GaussianSmootherComputationMode.EXPECTATION);
            case CANONICAL_ANALYTICAL:
                return createAnalyticalGradientEngine(latentProcess, process, observationModel, model, smootherMode);
            default:
                throw new IllegalArgumentException("Unsupported gradient mode: " + gradientMode);
        }
    }

    private static GradientEngine createAnalyticalGradientEngine(final LatentProcessModel latentProcess,
                                                                 final RepresentableProcess process,
                                                                 final GaussianObservationModel observationModel,
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
                "Expectation-form analytical gradients are not configured for process type: "
                        + process.getClass().getName());
    }

    private static GaussianSmootherResults createSmoother(final RepresentableProcess process,
                                                          final GaussianObservationModel observationModel,
                                                          final TimeSeriesModel model,
                                                          final GaussianSmootherComputationMode smootherMode) {
        final GaussianTransitionRepresentation transitionRepresentation =
                process.getRepresentation(GaussianTransitionRepresentation.class);
        switch (smootherMode) {
            case EXPECTATION:
                return new KalmanSmootherEngine(transitionRepresentation, observationModel, model.getTimeGrid());
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
                                                             final GaussianSmootherResults smoother,
                                                             final MatrixParameterInterface driftMatrix,
                                                             final dr.inference.model.MatrixParameterInterface diffusionMatrix,
                                                             final dr.inference.model.Parameter stationaryMean,
                                                             final MatrixParameter initialCovariance,
                                                             final int stateDimension) {
        final GradientFormulaBundle formulas = new GradientFormulaBundle(
                new SelectionMatrixGradientFormula(driftMatrix, stateDimension),
                processModel != null
                        ? new CanonicalSelectionMatrixGradientFormula(processModel, driftMatrix, stateDimension)
                        : new CanonicalSelectionMatrixGradientFormula(driftMatrix, stateDimension),
                new StationaryMeanGradientFormula(
                        processModel,
                        stationaryMean,
                        initialCovariance,
                        stateDimension),
                new CanonicalStationaryMeanGradientFormula(
                        processModel,
                        stationaryMean,
                        initialCovariance,
                        stateDimension),
                new DiffusionMatrixGradientFormula(
                        DiffusionMatrixParameterizationFactory.create(diffusionMatrix),
                        stateDimension),
                new CanonicalDiffusionMatrixGradientFormula(
                        processModel,
                        DiffusionMatrixParameterizationFactory.create(diffusionMatrix),
                        stateDimension));
        if (smoother instanceof CanonicalKalmanSmootherEngine) {
            return new CanonicalAnalyticalKalmanGradientEngine(
                    (CanonicalKalmanSmootherEngine) smoother,
                    formulas.canonicalSelection,
                    formulas.canonicalMean,
                    formulas.canonicalDiffusion);
        }
        return new AnalyticalKalmanGradientEngine(
                smoother,
                formulas.selection,
                formulas.mean,
                formulas.diffusion);
    }

    private static final class GradientFormulaBundle {
        final SelectionMatrixGradientFormula selection;
        final CanonicalSelectionMatrixGradientFormula canonicalSelection;
        final StationaryMeanGradientFormula mean;
        final CanonicalStationaryMeanGradientFormula canonicalMean;
        final DiffusionMatrixGradientFormula diffusion;
        final CanonicalDiffusionMatrixGradientFormula canonicalDiffusion;

        GradientFormulaBundle(final SelectionMatrixGradientFormula selection,
                              final CanonicalSelectionMatrixGradientFormula canonicalSelection,
                              final StationaryMeanGradientFormula mean,
                              final CanonicalStationaryMeanGradientFormula canonicalMean,
                              final DiffusionMatrixGradientFormula diffusion,
                              final CanonicalDiffusionMatrixGradientFormula canonicalDiffusion) {
            this.selection = selection;
            this.canonicalSelection = canonicalSelection;
            this.mean = mean;
            this.canonicalMean = canonicalMean;
            this.diffusion = diffusion;
            this.canonicalDiffusion = canonicalDiffusion;
        }
    }
}
