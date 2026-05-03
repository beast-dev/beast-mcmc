package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Exact block-space helper for Fréchet derivatives of exp(-D t).
 *
 * <p>The forward block integral is represented as a cached exact linear map on each
 * 1x1 / 2x2 block pair. This avoids numerical quadrature on the hot path and turns
 * repeated branch-local applications into tiny dense matrix-vector multiplies.</p>
 */
public final class BlockDiagonalFrechetHelper {

    private static final String QUADRATURE_DEBUG_PROPERTY =
            "beast.experimental.nativeUseQuadratureFrechet";
    private static final String EXACT_PLAN_DIAGNOSTICS_PROPERTY =
            "beast.debug.canonicalFrechetPlan";
    private static final String EXACT_PLAN_DIAGNOSTICS_REPORT_EVERY_PROPERTY =
            "beast.debug.canonicalFrechetPlanReportEvery";
    private static final boolean EXACT_PLAN_DIAGNOSTICS_ENABLED =
            Boolean.getBoolean(EXACT_PLAN_DIAGNOSTICS_PROPERTY);
    private static final int EXACT_PLAN_DIAGNOSTICS_REPORT_EVERY =
            Math.max(1, Integer.getInteger(EXACT_PLAN_DIAGNOSTICS_REPORT_EVERY_PROPERTY, 50));
    private static final AtomicLong EXACT_PLAN_CACHE_HITS = new AtomicLong();
    private static final AtomicLong EXACT_PLAN_CACHE_MISSES = new AtomicLong();
    private static final AtomicLong EXACT_PLAN_PARAMETER_MISSES = new AtomicLong();
    private static final AtomicLong EXACT_PLAN_TIME_MISSES = new AtomicLong();
    private static long exactPlanDiagnosticsCalls;
    private static long lastReportedParameterUpdates;
    private static long lastReportedTimeEvaluations;
    private static long lastReportedApplications;
    private static long lastReportedCacheHits;
    private static long lastReportedCacheMisses;
    private static long lastReportedParameterMisses;
    private static long lastReportedTimeMisses;
    private static long lastReportedEqualDiagonalEvaluations;
    private static long lastReportedGenericSolve4x4Calls;
    private static long lastReportedKernel1x1Evaluations;
    private static long lastReportedKernel1x2Evaluations;
    private static long lastReportedKernel2x1Evaluations;
    private static long lastReportedKernel2x2EqualDiagonalEvaluations;
    private static long lastReportedKernel2x2GenericEvaluations;
    private static long lastReportedCoefficientBivariateSmallRootEvaluations;
    private static long lastReportedCoefficientSmallLeftRootEvaluations;
    private static long lastReportedCoefficientSmallRightRootEvaluations;
    private static long lastReportedCoefficientDistinctEvaluations;
    private static long lastReportedCoefficientDistinctRealEvaluations;
    private static long lastReportedCoefficientDistinctLeftImagRightRealEvaluations;
    private static long lastReportedCoefficientDistinctLeftRealRightImagEvaluations;
    private static long lastReportedCoefficientDistinctBothImagEvaluations;

    private final int dim;
    private final BlockDiagonalExpSolver.BlockStructure structure;
    private final BlockDiagonalExpSolver expSolver;

    private final double[] blockDTransposeParams;
    private final double[] buf2a = new double[2];
    private final double[] buf4a = new double[4];
    private final double[] buf4b = new double[4];
    private final double[] buf4c = new double[4];
    private final DenseMatrix64F lastExp;
    private final BlockDiagonalFrechetExactPlan.Workspace workspace;

    private final double[] cachedPlanParams;
    private final BlockDiagonalFrechetExactPlan.Plan cachedPlan;
    private double cachedPlanT;
    private int cachedPlanParameterHash;
    private long cachedPlanParameterVersion;
    private boolean cachedPlanValid;

    public BlockDiagonalFrechetHelper(final BlockDiagonalExpSolver.BlockStructure structure) {
        this.structure = structure;
        this.dim = structure.getDim();
        this.expSolver = new BlockDiagonalExpSolver(structure);
        this.blockDTransposeParams = new double[3 * dim - 2];
        this.lastExp = new DenseMatrix64F(dim, dim);
        this.workspace = new BlockDiagonalFrechetExactPlan.Workspace();
        this.cachedPlanParams = new double[3 * dim - 2];
        this.cachedPlan = BlockDiagonalFrechetExactPlan.createPlan(structure);
        this.cachedPlanT = Double.NaN;
        this.cachedPlanParameterHash = 0;
        this.cachedPlanParameterVersion = 0L;
        this.cachedPlanValid = false;
    }

    public static void resetExactPlanInstrumentation() {
        BlockDiagonalFrechetExactPlan.resetInstrumentation();
        EXACT_PLAN_CACHE_HITS.set(0L);
        EXACT_PLAN_CACHE_MISSES.set(0L);
        EXACT_PLAN_PARAMETER_MISSES.set(0L);
        EXACT_PLAN_TIME_MISSES.set(0L);
        resetExactPlanDiagnosticsSnapshot();
    }

    public static long getExactPlanParameterUpdateCount() {
        return BlockDiagonalFrechetExactPlan.getParameterUpdateCount();
    }

    public static long getExactPlanTimeEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getTimeEvaluationCount();
    }

    public static long getExactPlanApplicationCount() {
        return BlockDiagonalFrechetExactPlan.getApplicationCount();
    }

    public static long getExactPlanCacheHitCount() {
        return EXACT_PLAN_CACHE_HITS.get();
    }

    public static long getExactPlanCacheMissCount() {
        return EXACT_PLAN_CACHE_MISSES.get();
    }

    public static long getExactPlanParameterMissCount() {
        return EXACT_PLAN_PARAMETER_MISSES.get();
    }

    public static long getExactPlanTimeMissCount() {
        return EXACT_PLAN_TIME_MISSES.get();
    }

    public static long getExactPlanEqualDiagonalEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getEqualDiagonalEvaluationCount();
    }

    public static long getExactPlanGenericSolve4x4CallCount() {
        return BlockDiagonalFrechetExactPlan.getGenericSolve4x4CallCount();
    }

    public static long getExactPlanKernel1x1EvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getKernel1x1EvaluationCount();
    }

    public static long getExactPlanKernel1x2EvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getKernel1x2EvaluationCount();
    }

    public static long getExactPlanKernel2x1EvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getKernel2x1EvaluationCount();
    }

    public static long getExactPlanKernel2x2EqualDiagonalEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getKernel2x2EqualDiagonalEvaluationCount();
    }

    public static long getExactPlanKernel2x2GenericEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getKernel2x2GenericEvaluationCount();
    }

    public static long getExactPlanCoefficientBivariateSmallRootEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientBivariateSmallRootEvaluationCount();
    }

    public static long getExactPlanCoefficientSmallLeftRootEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientSmallLeftRootEvaluationCount();
    }

    public static long getExactPlanCoefficientSmallRightRootEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientSmallRightRootEvaluationCount();
    }

    public static long getExactPlanCoefficientDistinctEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientDistinctEvaluationCount();
    }

    public static long getExactPlanCoefficientDistinctRealEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientDistinctRealEvaluationCount();
    }

    public static long getExactPlanCoefficientDistinctLeftImagRightRealEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientDistinctLeftImagRightRealEvaluationCount();
    }

    public static long getExactPlanCoefficientDistinctLeftRealRightImagEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientDistinctLeftRealRightImagEvaluationCount();
    }

    public static long getExactPlanCoefficientDistinctBothImagEvaluationCount() {
        return BlockDiagonalFrechetExactPlan.getCoefficientDistinctBothImagEvaluationCount();
    }

    public static void reportExactPlanDiagnosticsIfEnabled() {
        if (!EXACT_PLAN_DIAGNOSTICS_ENABLED) {
            return;
        }
        exactPlanDiagnosticsCalls++;
        if (exactPlanDiagnosticsCalls % EXACT_PLAN_DIAGNOSTICS_REPORT_EVERY != 0L) {
            return;
        }
        reportExactPlanDiagnosticsAndResetSnapshot();
    }

    public DenseMatrix64F computeExpInDBasis(final double[] blockDParams,
                                             final double t) {
        expSolver.compute(blockDParams, t, lastExp);
        return lastExp;
    }

    public void frechetAdjointExpInDBasis(final double[] blockDParams,
                                          final DenseMatrix64F xD,
                                          final double t,
                                          final DenseMatrix64F out) {
        validateInputs(blockDParams, xD, out);
        buildTransposeParams(blockDParams, blockDTransposeParams);
        computeForwardFrechetInDBasis(blockDTransposeParams, xD, t, out);
        CommonOps.transpose(out);
        CommonOps.scale(-t, out);
    }

    public void computeForwardFrechetInDBasis(final double[] blockDParams,
                                              final DenseMatrix64F eD,
                                              final double t,
                                              final DenseMatrix64F out) {
        validateInputs(blockDParams, eD, out);
        Arrays.fill(out.data, 0.0);
        if (Boolean.getBoolean(QUADRATURE_DEBUG_PROPERTY)) {
            computeForwardFrechetInDBasisQuadrature(blockDParams, eD, t, out);
            return;
        }
        getOrBuildPlan(blockDParams, t).apply(eD, out, workspace);
    }

    private void buildTransposeParams(final double[] src, final double[] dst) {
        final int diagLen = dim;
        final int offLen = dim - 1;
        System.arraycopy(src, 0, dst, 0, diagLen);
        System.arraycopy(src, diagLen + offLen, dst, diagLen, offLen);
        System.arraycopy(src, diagLen, dst, diagLen + offLen, offLen);
    }

    private void validateInputs(final double[] blockDParams,
                                final DenseMatrix64F in,
                                final DenseMatrix64F out) {
        BlockDiagonalExpSolver.validateCompressedParams(blockDParams, dim);
        validateMatrix(in, "input");
        validateMatrix(out, "out");
    }

    private void validateMatrix(final DenseMatrix64F x, final String name) {
        if (x.numRows != dim || x.numCols != dim) {
            throw new IllegalArgumentException(
                    name + " must be " + dim + "x" + dim + " but is " + x.numRows + "x" + x.numCols);
        }
    }

    private void computeForwardFrechetInDBasisQuadrature(final double[] blockDParams,
                                                         final DenseMatrix64F eD,
                                                         final double t,
                                                         final DenseMatrix64F out) {
        for (int bi = 0; bi < structure.getNumBlocks(); bi++) {
            final int startI = structure.getBlockStart(bi);
            final int sizeI = structure.getBlockSize(bi);
            for (int bj = 0; bj < structure.getNumBlocks(); bj++) {
                final int startJ = structure.getBlockStart(bj);
                final int sizeJ = structure.getBlockSize(bj);
                computeBlockFrechetIntegralQuadrature(blockDParams, eD, t, startI, sizeI, startJ, sizeJ, out);
            }
        }
    }

    private void computeBlockFrechetIntegralQuadrature(final double[] blockDParams,
                                                       final DenseMatrix64F eTilde,
                                                       final double t,
                                                       final int startI,
                                                       final int sizeI,
                                                       final int startJ,
                                                       final int sizeJ,
                                                       final DenseMatrix64F out) {
        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        if (sizeI == 1 && sizeJ == 1) {
            final double lambdaI = -t * blockDParams[startI];
            final double lambdaJ = -t * blockDParams[startJ];
            final double eIJ = eTilde.get(startI, startJ);
            final double value;
            if (Math.abs(lambdaI - lambdaJ) < 1.0e-10) {
                value = eIJ * Math.exp(lambdaI);
            } else {
                value = eIJ * (Math.exp(lambdaI) - Math.exp(lambdaJ)) / (lambdaI - lambdaJ);
            }
            out.set(startI, startJ, value);
            return;
        }

        if (sizeI == 1 && sizeJ == 2) {
            final double lambdaI = -t * blockDParams[startI];
            final double a = -t * blockDParams[startJ];
            final double d = -t * blockDParams[startJ + 1];
            final double b = -t * blockDParams[upperOffset + startJ];
            final double c = -t * blockDParams[lowerOffset + startJ];
            final double e1 = eTilde.get(startI, startJ);
            final double e2 = eTilde.get(startI, startJ + 1);

            BlockDiagonalFrechetIntegrator.integrate1x1_2x2(
                    lambdaI, a, b, c, d, e1, e2, buf2a, buf4a);
            out.set(startI, startJ, buf2a[0]);
            out.set(startI, startJ + 1, buf2a[1]);
            return;
        }

        if (sizeI == 2 && sizeJ == 1) {
            final double a = -t * blockDParams[startI];
            final double d = -t * blockDParams[startI + 1];
            final double b = -t * blockDParams[upperOffset + startI];
            final double c = -t * blockDParams[lowerOffset + startI];
            final double lambdaJ = -t * blockDParams[startJ];
            final double e1 = eTilde.get(startI, startJ);
            final double e2 = eTilde.get(startI + 1, startJ);

            BlockDiagonalFrechetIntegrator.integrate2x2_1x1(
                    a, b, c, d, lambdaJ, e1, e2, buf2a, buf4a);
            out.set(startI, startJ, buf2a[0]);
            out.set(startI + 1, startJ, buf2a[1]);
            return;
        }

        final double ai = -t * blockDParams[startI];
        final double di = -t * blockDParams[startI + 1];
        final double bi = -t * blockDParams[upperOffset + startI];
        final double ci = -t * blockDParams[lowerOffset + startI];

        final double aj = -t * blockDParams[startJ];
        final double dj = -t * blockDParams[startJ + 1];
        final double bj = -t * blockDParams[upperOffset + startJ];
        final double cj = -t * blockDParams[lowerOffset + startJ];

        final double e11 = eTilde.get(startI, startJ);
        final double e12 = eTilde.get(startI, startJ + 1);
        final double e21 = eTilde.get(startI + 1, startJ);
        final double e22 = eTilde.get(startI + 1, startJ + 1);

        BlockDiagonalFrechetIntegrator.integrate2x2_2x2(
                ai, bi, ci, di,
                aj, bj, cj, dj,
                e11, e12, e21, e22,
                buf4a, buf4b, buf4c
        );
        out.set(startI, startJ, buf4a[0]);
        out.set(startI, startJ + 1, buf4a[1]);
        out.set(startI + 1, startJ, buf4a[2]);
        out.set(startI + 1, startJ + 1, buf4a[3]);
    }

    private BlockDiagonalFrechetExactPlan.Plan getOrBuildPlan(final double[] blockDParams,
                                                              final double t) {
        final int parameterHash = Arrays.hashCode(blockDParams);
        final boolean parametersChanged = !cachedPlanValid
                || parameterHash != cachedPlanParameterHash
                || !Arrays.equals(blockDParams, cachedPlanParams);
        final boolean timeChanged = !cachedPlanValid
                || Double.doubleToLongBits(t) != Double.doubleToLongBits(cachedPlanT);
        if (!parametersChanged && !timeChanged) {
            EXACT_PLAN_CACHE_HITS.incrementAndGet();
            return cachedPlan;
        }

        EXACT_PLAN_CACHE_MISSES.incrementAndGet();
        if (parametersChanged) {
            EXACT_PLAN_PARAMETER_MISSES.incrementAndGet();
        }
        if (timeChanged) {
            EXACT_PLAN_TIME_MISSES.incrementAndGet();
        }

        if (parametersChanged || timeChanged) {
            if (parametersChanged) {
                System.arraycopy(blockDParams, 0, cachedPlanParams, 0, blockDParams.length);
                cachedPlanParameterHash = parameterHash;
                cachedPlanParameterVersion++;
                cachedPlan.updateParameters(blockDParams);
            }
            cachedPlan.evaluate(t);
            cachedPlanT = t;
            cachedPlanValid = true;
        }
        return cachedPlan;
    }

    private static void reportExactPlanDiagnosticsAndResetSnapshot() {
        final long parameterUpdates = getExactPlanParameterUpdateCount();
        final long timeEvaluations = getExactPlanTimeEvaluationCount();
        final long applications = getExactPlanApplicationCount();
        final long cacheHits = getExactPlanCacheHitCount();
        final long cacheMisses = getExactPlanCacheMissCount();
        final long parameterMisses = getExactPlanParameterMissCount();
        final long timeMisses = getExactPlanTimeMissCount();
        final long equalDiagonalEvaluations = getExactPlanEqualDiagonalEvaluationCount();
        final long genericSolve4x4Calls = getExactPlanGenericSolve4x4CallCount();
        final long kernel1x1Evaluations = getExactPlanKernel1x1EvaluationCount();
        final long kernel1x2Evaluations = getExactPlanKernel1x2EvaluationCount();
        final long kernel2x1Evaluations = getExactPlanKernel2x1EvaluationCount();
        final long kernel2x2EqualDiagonalEvaluations = getExactPlanKernel2x2EqualDiagonalEvaluationCount();
        final long kernel2x2GenericEvaluations = getExactPlanKernel2x2GenericEvaluationCount();
        final long coefficientBivariateSmallRootEvaluations =
                getExactPlanCoefficientBivariateSmallRootEvaluationCount();
        final long coefficientSmallLeftRootEvaluations = getExactPlanCoefficientSmallLeftRootEvaluationCount();
        final long coefficientSmallRightRootEvaluations = getExactPlanCoefficientSmallRightRootEvaluationCount();
        final long coefficientDistinctEvaluations = getExactPlanCoefficientDistinctEvaluationCount();
        final long coefficientDistinctRealEvaluations = getExactPlanCoefficientDistinctRealEvaluationCount();
        final long coefficientDistinctLeftImagRightRealEvaluations =
                getExactPlanCoefficientDistinctLeftImagRightRealEvaluationCount();
        final long coefficientDistinctLeftRealRightImagEvaluations =
                getExactPlanCoefficientDistinctLeftRealRightImagEvaluationCount();
        final long coefficientDistinctBothImagEvaluations =
                getExactPlanCoefficientDistinctBothImagEvaluationCount();

        final long deltaParameterUpdates = parameterUpdates - lastReportedParameterUpdates;
        final long deltaTimeEvaluations = timeEvaluations - lastReportedTimeEvaluations;
        final long deltaApplications = applications - lastReportedApplications;
        final long deltaCacheHits = cacheHits - lastReportedCacheHits;
        final long deltaCacheMisses = cacheMisses - lastReportedCacheMisses;
        final long deltaParameterMisses = parameterMisses - lastReportedParameterMisses;
        final long deltaTimeMisses = timeMisses - lastReportedTimeMisses;
        final long deltaEqualDiagonalEvaluations =
                equalDiagonalEvaluations - lastReportedEqualDiagonalEvaluations;
        final long deltaGenericSolve4x4Calls = genericSolve4x4Calls - lastReportedGenericSolve4x4Calls;
        final long deltaKernel1x1Evaluations = kernel1x1Evaluations - lastReportedKernel1x1Evaluations;
        final long deltaKernel1x2Evaluations = kernel1x2Evaluations - lastReportedKernel1x2Evaluations;
        final long deltaKernel2x1Evaluations = kernel2x1Evaluations - lastReportedKernel2x1Evaluations;
        final long deltaKernel2x2EqualDiagonalEvaluations =
                kernel2x2EqualDiagonalEvaluations - lastReportedKernel2x2EqualDiagonalEvaluations;
        final long deltaKernel2x2GenericEvaluations =
                kernel2x2GenericEvaluations - lastReportedKernel2x2GenericEvaluations;
        final long deltaCoefficientBivariateSmallRootEvaluations =
                coefficientBivariateSmallRootEvaluations - lastReportedCoefficientBivariateSmallRootEvaluations;
        final long deltaCoefficientSmallLeftRootEvaluations =
                coefficientSmallLeftRootEvaluations - lastReportedCoefficientSmallLeftRootEvaluations;
        final long deltaCoefficientSmallRightRootEvaluations =
                coefficientSmallRightRootEvaluations - lastReportedCoefficientSmallRightRootEvaluations;
        final long deltaCoefficientDistinctEvaluations =
                coefficientDistinctEvaluations - lastReportedCoefficientDistinctEvaluations;
        final long deltaCoefficientDistinctRealEvaluations =
                coefficientDistinctRealEvaluations - lastReportedCoefficientDistinctRealEvaluations;
        final long deltaCoefficientDistinctLeftImagRightRealEvaluations =
                coefficientDistinctLeftImagRightRealEvaluations
                        - lastReportedCoefficientDistinctLeftImagRightRealEvaluations;
        final long deltaCoefficientDistinctLeftRealRightImagEvaluations =
                coefficientDistinctLeftRealRightImagEvaluations
                        - lastReportedCoefficientDistinctLeftRealRightImagEvaluations;
        final long deltaCoefficientDistinctBothImagEvaluations =
                coefficientDistinctBothImagEvaluations
                        - lastReportedCoefficientDistinctBothImagEvaluations;
        final long cacheRequests = deltaCacheHits + deltaCacheMisses;

        final String message = String.format(Locale.US,
                "[canonical-frechet-plan] gradientReports=%d window.applications=%d "
                        + "window.timeEvaluations=%d window.parameterUpdates=%d "
                        + "window.cache=requests/hits/misses=%d/%d/%d hitRate=%.1f%% "
                        + "window.missCauses=parameter/time=%d/%d "
                        + "window.kernels=1x1/1x2/2x1/2x2equal/2x2generic=%d/%d/%d/%d/%d "
                        + "window.equalDiagCoeff=bivarSmall/smallLeft/smallRight/distinct=%d/%d/%d/%d "
                        + "window.distinctShapes=real/leftImagRightReal/leftRealRightImag/bothImag=%d/%d/%d/%d "
                        + "window.equalDiagonal=%d window.genericSolve4x4=%d",
                exactPlanDiagnosticsCalls,
                deltaApplications,
                deltaTimeEvaluations,
                deltaParameterUpdates,
                cacheRequests,
                deltaCacheHits,
                deltaCacheMisses,
                cacheRequests == 0L ? 0.0 : 100.0 * deltaCacheHits / cacheRequests,
                deltaParameterMisses,
                deltaTimeMisses,
                deltaKernel1x1Evaluations,
                deltaKernel1x2Evaluations,
                deltaKernel2x1Evaluations,
                deltaKernel2x2EqualDiagonalEvaluations,
                deltaKernel2x2GenericEvaluations,
                deltaCoefficientBivariateSmallRootEvaluations,
                deltaCoefficientSmallLeftRootEvaluations,
                deltaCoefficientSmallRightRootEvaluations,
                deltaCoefficientDistinctEvaluations,
                deltaCoefficientDistinctRealEvaluations,
                deltaCoefficientDistinctLeftImagRightRealEvaluations,
                deltaCoefficientDistinctLeftRealRightImagEvaluations,
                deltaCoefficientDistinctBothImagEvaluations,
                deltaEqualDiagonalEvaluations,
                deltaGenericSolve4x4Calls);
        System.err.println(message);

        lastReportedParameterUpdates = parameterUpdates;
        lastReportedTimeEvaluations = timeEvaluations;
        lastReportedApplications = applications;
        lastReportedCacheHits = cacheHits;
        lastReportedCacheMisses = cacheMisses;
        lastReportedParameterMisses = parameterMisses;
        lastReportedTimeMisses = timeMisses;
        lastReportedEqualDiagonalEvaluations = equalDiagonalEvaluations;
        lastReportedGenericSolve4x4Calls = genericSolve4x4Calls;
        lastReportedKernel1x1Evaluations = kernel1x1Evaluations;
        lastReportedKernel1x2Evaluations = kernel1x2Evaluations;
        lastReportedKernel2x1Evaluations = kernel2x1Evaluations;
        lastReportedKernel2x2EqualDiagonalEvaluations = kernel2x2EqualDiagonalEvaluations;
        lastReportedKernel2x2GenericEvaluations = kernel2x2GenericEvaluations;
        lastReportedCoefficientBivariateSmallRootEvaluations = coefficientBivariateSmallRootEvaluations;
        lastReportedCoefficientSmallLeftRootEvaluations = coefficientSmallLeftRootEvaluations;
        lastReportedCoefficientSmallRightRootEvaluations = coefficientSmallRightRootEvaluations;
        lastReportedCoefficientDistinctEvaluations = coefficientDistinctEvaluations;
        lastReportedCoefficientDistinctRealEvaluations = coefficientDistinctRealEvaluations;
        lastReportedCoefficientDistinctLeftImagRightRealEvaluations =
                coefficientDistinctLeftImagRightRealEvaluations;
        lastReportedCoefficientDistinctLeftRealRightImagEvaluations =
                coefficientDistinctLeftRealRightImagEvaluations;
        lastReportedCoefficientDistinctBothImagEvaluations = coefficientDistinctBothImagEvaluations;
    }

    private static void resetExactPlanDiagnosticsSnapshot() {
        exactPlanDiagnosticsCalls = 0L;
        lastReportedParameterUpdates = 0L;
        lastReportedTimeEvaluations = 0L;
        lastReportedApplications = 0L;
        lastReportedCacheHits = 0L;
        lastReportedCacheMisses = 0L;
        lastReportedParameterMisses = 0L;
        lastReportedTimeMisses = 0L;
        lastReportedEqualDiagonalEvaluations = 0L;
        lastReportedGenericSolve4x4Calls = 0L;
        lastReportedKernel1x1Evaluations = 0L;
        lastReportedKernel1x2Evaluations = 0L;
        lastReportedKernel2x1Evaluations = 0L;
        lastReportedKernel2x2EqualDiagonalEvaluations = 0L;
        lastReportedKernel2x2GenericEvaluations = 0L;
        lastReportedCoefficientBivariateSmallRootEvaluations = 0L;
        lastReportedCoefficientSmallLeftRootEvaluations = 0L;
        lastReportedCoefficientSmallRightRootEvaluations = 0L;
        lastReportedCoefficientDistinctEvaluations = 0L;
        lastReportedCoefficientDistinctRealEvaluations = 0L;
        lastReportedCoefficientDistinctLeftImagRightRealEvaluations = 0L;
        lastReportedCoefficientDistinctLeftRealRightImagEvaluations = 0L;
        lastReportedCoefficientDistinctBothImagEvaluations = 0L;
    }
}
