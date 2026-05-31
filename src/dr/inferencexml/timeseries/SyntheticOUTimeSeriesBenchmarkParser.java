package dr.inferencexml.timeseries;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.runtime.CompositeGradientWrtParameterProvider;
import dr.inference.timeseries.runtime.TimeSeriesBenchmark;
import dr.inference.timeseries.runtime.TimeSeriesGradient;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * XML parser for the deterministic OU time-series benchmark used for Stan comparisons.
 */
public class SyntheticOUTimeSeriesBenchmarkParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "ouTimeSeriesBenchmark";
    public static final String MODE = "mode";
    public static final String STATE_DIMENSION = "stateDimension";
    public static final String TIME_COUNT = "timeCount";
    public static final String TIME_STEP = "timeStep";
    public static final String SIGMA_OBS = "sigmaObs";
    public static final String GRID = "grid";
    public static final String WARMUP = "warmup";
    public static final String N_WARMUP = "nWarmup";
    public static final String TIMED_ITERATIONS = "timedIterations";
    public static final String N_TIMED = "nTimed";
    public static final String ITERATION_COUNT = "iterationCount";
    public static final String PERTURB_EACH_ITERATION = "perturbEachIteration";
    public static final String PERTURBATION = "perturbation";
    public static final String PERTURB_INDEX = "perturbIndex";
    public static final String SELECTION_GRADIENT = "selectionGradient";
    public static final String DIFFUSION_GRADIENT = "diffusionGradient";
    public static final String MEAN_GRADIENT = "meanGradient";

    private static final String GRID_UNIFORM = "uniform";
    private static final String GRID_IRREGULAR = "irregular";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        final TimeSeriesBenchmark.Mode mode;
        try {
            mode = TimeSeriesBenchmark.Mode.parse(xo.getAttribute(MODE, "logLikelihood"));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }

        final int stateDimension = xo.getAttribute(STATE_DIMENSION, 10);
        final int timeCount = xo.getAttribute(TIME_COUNT, 300);
        final double timeStep = xo.getAttribute(TIME_STEP, 0.05);
        final double sigmaObs = xo.getAttribute(SIGMA_OBS, 0.10);
        final String gridType = xo.getAttribute(GRID, GRID_IRREGULAR);
        final int warmup = getAliasedIntegerAttribute(xo, WARMUP, N_WARMUP, 20);
        final int timedIterations = getTimedIterations(xo);
        final boolean perturbEachIteration = xo.getAttribute(PERTURB_EACH_ITERATION, true);
        final double perturbation = xo.getAttribute(PERTURBATION, 1.0e-9);
        final int perturbIndex = xo.getAttribute(PERTURB_INDEX, 0);
        final boolean selectionGradient = xo.getAttribute(SELECTION_GRADIENT, true);
        final boolean diffusionGradient = xo.getAttribute(DIFFUSION_GRADIENT, false);
        final boolean meanGradient = xo.getAttribute(MEAN_GRADIENT, false);

        validateSettings(stateDimension, timeCount, timeStep, sigmaObs, gridType);
        validateGradientSettings(mode, selectionGradient, diffusionGradient, meanGradient);

        final BenchmarkModel benchmarkModel = makeBenchmarkModel(
                id,
                stateDimension,
                timeCount,
                timeStep,
                sigmaObs,
                gridType,
                selectionGradient,
                diffusionGradient,
                meanGradient);
        final String description = String.format(Locale.US,
                "K=%d  T=%d  grid=%s  sigmaObs=%.4f  uniqueDt=%d  gradientComponents=%s",
                stateDimension, timeCount, gridType, sigmaObs, benchmarkModel.uniqueDeltaCount,
                benchmarkModel.gradientComponents);

        return new TimeSeriesBenchmark(
                id,
                mode,
                benchmarkModel.likelihood,
                benchmarkModel.gradient,
                warmup,
                timedIterations,
                perturbEachIteration,
                benchmarkModel.perturbParameter,
                perturbIndex,
                perturbation,
                description);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newStringRule(MODE, true,
                    "Benchmark mode: logLikelihood/logL or gradient/grad. Defaults to logLikelihood."),
            AttributeRule.newIntegerRule(STATE_DIMENSION, true),
            AttributeRule.newIntegerRule(TIME_COUNT, true),
            AttributeRule.newDoubleRule(TIME_STEP, true),
            AttributeRule.newDoubleRule(SIGMA_OBS, true),
            new StringAttributeRule(GRID,
                    "Time grid: irregular or uniform. Defaults to irregular.",
                    new String[] { GRID_IRREGULAR, GRID_UNIFORM }, true),
            AttributeRule.newIntegerRule(WARMUP, true),
            AttributeRule.newIntegerRule(N_WARMUP, true),
            AttributeRule.newIntegerRule(TIMED_ITERATIONS, true),
            AttributeRule.newIntegerRule(N_TIMED, true),
            AttributeRule.newIntegerRule(ITERATION_COUNT, true),
            AttributeRule.newBooleanRule(PERTURB_EACH_ITERATION, true),
            AttributeRule.newDoubleRule(PERTURBATION, true),
            AttributeRule.newIntegerRule(PERTURB_INDEX, true),
            AttributeRule.newBooleanRule(SELECTION_GRADIENT, true),
            AttributeRule.newBooleanRule(DIFFUSION_GRADIENT, true),
            AttributeRule.newBooleanRule(MEAN_GRADIENT, true)
    };

    @Override
    public String getParserDescription() {
        return "Runs a deterministic synthetic OU time-series benchmark for log-likelihood or gradient evaluation.";
    }

    @Override
    public Class getReturnType() {
        return TimeSeriesBenchmark.class;
    }

    private static int getTimedIterations(final XMLObject xo) throws XMLParseException {
        if (xo.hasAttribute(TIMED_ITERATIONS)) {
            return xo.getIntegerAttribute(TIMED_ITERATIONS);
        }
        if (xo.hasAttribute(N_TIMED)) {
            return xo.getIntegerAttribute(N_TIMED);
        }
        if (xo.hasAttribute(ITERATION_COUNT)) {
            return xo.getIntegerAttribute(ITERATION_COUNT);
        }
        return 100;
    }

    private static int getAliasedIntegerAttribute(final XMLObject xo,
                                                  final String primary,
                                                  final String alias,
                                                  final int defaultValue) throws XMLParseException {
        if (xo.hasAttribute(primary)) {
            return xo.getIntegerAttribute(primary);
        }
        if (xo.hasAttribute(alias)) {
            return xo.getIntegerAttribute(alias);
        }
        return defaultValue;
    }

    private static void validateSettings(final int stateDimension,
                                         final int timeCount,
                                         final double timeStep,
                                         final double sigmaObs,
                                         final String gridType) throws XMLParseException {
        if (stateDimension < 2 || (stateDimension & 1) != 0) {
            throw new XMLParseException(STATE_DIMENSION + " must be an even integer >= 2");
        }
        if (timeCount < 2) {
            throw new XMLParseException(TIME_COUNT + " must be at least 2");
        }
        if (!(timeStep > 0.0)) {
            throw new XMLParseException(TIME_STEP + " must be positive");
        }
        if (!(sigmaObs > 0.0)) {
            throw new XMLParseException(SIGMA_OBS + " must be positive");
        }
        if (!GRID_IRREGULAR.equalsIgnoreCase(gridType)
                && !GRID_UNIFORM.equalsIgnoreCase(gridType)) {
            throw new XMLParseException(GRID + " must be '" + GRID_IRREGULAR
                    + "' or '" + GRID_UNIFORM + "'");
        }
    }

    private static void validateGradientSettings(final TimeSeriesBenchmark.Mode mode,
                                                 final boolean selectionGradient,
                                                 final boolean diffusionGradient,
                                                 final boolean meanGradient) throws XMLParseException {
        if (mode == TimeSeriesBenchmark.Mode.GRADIENT
                && !selectionGradient
                && !diffusionGradient
                && !meanGradient) {
            throw new XMLParseException("At least one of " + SELECTION_GRADIENT + ", "
                    + DIFFUSION_GRADIENT + ", or " + MEAN_GRADIENT
                    + " must be true in gradient mode");
        }
    }

    private static BenchmarkModel makeBenchmarkModel(final String id,
                                                     final int stateDimension,
                                                     final int timeCount,
                                                     final double timeStep,
                                                     final double sigmaObs,
                                                     final String gridType,
                                                     final boolean selectionGradient,
                                                     final boolean diffusionGradient,
                                                     final boolean meanGradient) {
        final OrthogonalBlockDiagonalPolarStableMatrixParameter drift =
                makeOrthogonalBlockDrift(id + ".A", stateDimension);
        final OUProcessModel process = new OUProcessModel(
                id + ".process",
                stateDimension,
                drift,
                makeSpdMatrix(id + ".Q", stateDimension, 1.0, 0.12),
                makeMean(id + ".mu", stateDimension),
                makeSpdMatrix(id + ".P0", stateDimension, 1.0, 0.12));
        final TimeGrid grid = GRID_UNIFORM.equalsIgnoreCase(gridType)
                ? new UniformTimeGrid(timeCount, 0.0, timeStep)
                : makeIrregularGrid(timeCount, timeStep);
        final GaussianObservationModel observation =
                makeObservation(id + ".obs", stateDimension, timeCount, sigmaObs);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(
                id + ".model",
                new OUTimeSeriesProcessAdapter(process),
                observation,
                grid);
        final TimeSeriesLikelihood likelihood = GaussianTimeSeriesLikelihoodFactory.create(
                id + ".likelihood",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
        final GradientWrtParameterProvider gradient = makeGradientProvider(
                id,
                likelihood,
                drift,
                process,
                selectionGradient,
                diffusionGradient,
                meanGradient);

        return new BenchmarkModel(
                likelihood,
                gradient,
                drift.getRhoParameter(),
                countUniqueDeltas(grid, timeCount),
                gradientComponentDescription(
                        selectionGradient,
                        diffusionGradient,
                        meanGradient));
    }

    private static GradientWrtParameterProvider makeGradientProvider(
            final String id,
            final TimeSeriesLikelihood likelihood,
            final OrthogonalBlockDiagonalPolarStableMatrixParameter drift,
            final OUProcessModel process,
            final boolean selectionGradient,
            final boolean diffusionGradient,
            final boolean meanGradient) {
        final List<GradientWrtParameterProvider> providers =
                new ArrayList<GradientWrtParameterProvider>();
        if (selectionGradient) {
            providers.add(new TimeSeriesGradient(likelihood, drift.getParameter()));
        }
        if (diffusionGradient) {
            providers.add(new TimeSeriesGradient(likelihood, process.getDiffusionMatrix()));
        }
        if (meanGradient) {
            providers.add(new TimeSeriesGradient(likelihood, process.getStationaryMeanParameter()));
        }
        if (providers.isEmpty()) {
            return null;
        }
        if (providers.size() == 1) {
            return providers.get(0);
        }
        return new CompositeGradientWrtParameterProvider(id + ".fullGradient", providers);
    }

    private static String gradientComponentDescription(final boolean selectionGradient,
                                                       final boolean diffusionGradient,
                                                       final boolean meanGradient) {
        final StringBuilder builder = new StringBuilder();
        appendGradientComponent(builder, selectionGradient, "selection");
        appendGradientComponent(builder, diffusionGradient, "diffusion");
        appendGradientComponent(builder, meanGradient, "mean");
        return builder.length() == 0 ? "none" : builder.toString();
    }

    private static void appendGradientComponent(final StringBuilder builder,
                                                final boolean enabled,
                                                final String label) {
        if (enabled) {
            if (builder.length() > 0) {
                builder.append(',');
            }
            builder.append(label);
        }
    }

    private static OrthogonalBlockDiagonalPolarStableMatrixParameter makeOrthogonalBlockDrift(
            final String name,
            final int stateDimension) {
        final int nAngles = stateDimension * (stateDimension - 1) / 2;
        final int nBlocks = stateDimension / 2;

        final double[] anglesData = new double[nAngles];
        for (int i = 0; i < nAngles; ++i) {
            anglesData[i] = 0.15 * Math.sin(0.7 * i + 0.3);
        }
        final Parameter angles = new Parameter.Default(name + ".angles", anglesData);
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter(name + ".rotation", stateDimension, angles);

        final double[] rhoData = new double[nBlocks];
        final double[] thetaData = new double[nBlocks];
        final double[] tData = new double[nBlocks];
        for (int b = 0; b < nBlocks; ++b) {
            rhoData[b] = 0.60 + 0.15 * Math.sin(1.1 * b);
            thetaData[b] = 0.25 + 0.10 * Math.cos(0.9 * b + 0.5);
            tData[b] = 0.05 * Math.sin(1.3 * b + 1.0);
        }

        return new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                name,
                rotation,
                new Parameter.Default(0),
                new Parameter.Default(name + ".rho", rhoData),
                new Parameter.Default(name + ".theta", thetaData),
                new Parameter.Default(name + ".t", tData));
    }

    private static MatrixParameter makeSpdMatrix(final String name,
                                                 final int dimension,
                                                 final double diagVal,
                                                 final double offScale) {
        final double[][] data = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            data[i][i] = diagVal;
            for (int j = i + 1; j < dimension; ++j) {
                final double v = offScale * Math.sin(i + j + 1.0);
                data[i][j] = v;
                data[j][i] = v;
            }
        }
        final double loading = requiredDiagonalLoading(data);
        if (loading > 0.0) {
            for (int i = 0; i < dimension; ++i) {
                data[i][i] += loading;
            }
        }
        return makeMatrix(name, data);
    }

    private static double requiredDiagonalLoading(final double[][] matrix) {
        double loading = 0.0;
        for (int attempt = 0; attempt < 14; ++attempt) {
            if (isPositiveDefiniteWithLoading(matrix, loading)) {
                return loading;
            }
            loading = loading == 0.0 ? 1.0e-10 : loading * 10.0;
        }

        double minGershgorinMargin = Double.POSITIVE_INFINITY;
        for (int i = 0; i < matrix.length; ++i) {
            double radius = 0.0;
            for (int j = 0; j < matrix.length; ++j) {
                if (i != j) {
                    radius += Math.abs(matrix[i][j]);
                }
            }
            minGershgorinMargin = Math.min(minGershgorinMargin, matrix[i][i] - radius);
        }
        return minGershgorinMargin > 1.0e-8 ? 0.0 : 1.0e-8 - minGershgorinMargin;
    }

    private static boolean isPositiveDefiniteWithLoading(final double[][] matrix,
                                                         final double loading) {
        final int dimension = matrix.length;
        final double[][] cholesky = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i][j];
                if (i == j) {
                    sum += loading;
                }
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 1.0e-12)) {
                        return false;
                    }
                    cholesky[i][j] = Math.sqrt(sum);
                } else {
                    cholesky[i][j] = sum / cholesky[j][j];
                }
            }
        }
        return true;
    }

    private static Parameter makeMean(final String name, final int dimension) {
        final double[] data = new double[dimension];
        for (int i = 0; i < dimension; ++i) {
            data[i] = 0.10 * Math.sin(0.5 * i + 0.2);
        }
        return new Parameter.Default(name, data);
    }

    private static GaussianObservationModel makeObservation(final String name,
                                                           final int dimension,
                                                           final int timeCount,
                                                           final double sigmaObs) {
        final double[][] hData = new double[dimension][dimension];
        final double[][] rData = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            hData[i][i] = 1.0;
            rData[i][i] = sigmaObs * sigmaObs;
        }

        final double[][] y = new double[dimension][timeCount];
        for (int k = 0; k < dimension; ++k) {
            for (int t = 0; t < timeCount; ++t) {
                y[k][t] = 0.25 * Math.sin(0.08 * t + 0.5 * k)
                        + 0.06 * Math.cos(0.21 * t - 0.3 * k);
            }
        }

        return new GaussianObservationModel(
                name,
                dimension,
                makeMatrix(name + ".H", hData),
                makeMatrix(name + ".R", rData),
                makeMatrix(name + ".Y", y));
    }

    private static MatrixParameter makeMatrix(final String name, final double[][] values) {
        final MatrixParameter matrix = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < values[0].length; ++j) {
                matrix.setParameterValue(i, j, values[i][j]);
            }
        }
        return matrix;
    }

    private static IrregularTimeGrid makeIrregularGrid(final int timeCount, final double timeStep) {
        final double[] times = new double[timeCount];
        times[0] = 0.0;
        for (int i = 1; i < timeCount; ++i) {
            final double jitter = 0.4 * timeStep * Math.sin(2.7183 * i + 1.4142);
            times[i] = times[i - 1] + timeStep + jitter;
        }
        return new IrregularTimeGrid(times);
    }

    private static int countUniqueDeltas(final TimeGrid grid, final int timeCount) {
        final double[] unique = new double[timeCount - 1];
        int uniqueCount = 0;
        for (int i = 1; i < timeCount; ++i) {
            final double delta = grid.getTime(i) - grid.getTime(i - 1);
            boolean found = false;
            for (int j = 0; j < uniqueCount; ++j) {
                if (Math.abs(delta - unique[j]) < 1.0e-14) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                unique[uniqueCount++] = delta;
            }
        }
        return uniqueCount;
    }

    private static final class BenchmarkModel {
        final TimeSeriesLikelihood likelihood;
        final GradientWrtParameterProvider gradient;
        final Parameter perturbParameter;
        final int uniqueDeltaCount;
        final String gradientComponents;

        BenchmarkModel(final TimeSeriesLikelihood likelihood,
                       final GradientWrtParameterProvider gradient,
                       final Parameter perturbParameter,
                       final int uniqueDeltaCount,
                       final String gradientComponents) {
            this.likelihood = likelihood;
            this.gradient = gradient;
            this.perturbParameter = perturbParameter;
            this.uniqueDeltaCount = uniqueDeltaCount;
            this.gradientComponents = gradientComponents;
        }
    }
}
