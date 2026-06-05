/*
 * SericolaSeriesMarkovRewardFastModel.java
 *
 * Continuous reward density ONLY, in proportion space:
 *   alpha_i in [0,1] (reward-rate proportions)
 *   rewardProportion in (0,1)
 *
 * Physical mapping (handled outside this class):
 *   a_i = L + w * alpha_i
 *   r   = L + w * rewardProportion
 *
 * Assumptions:
 *  - Assumes alpha values are (after sorting) nondecreasing, and strictly increasing
 *    where intervals are used (ties are not allowed where h requires 1/(alpha_h-alpha_{h-1})).
 *  - Endpoints alpha_{i_L}=0 and alpha_{i_U}=1 are allowed.
 *  - rewardProportion is expected in [min alpha, max alpha]
 *  - This returns ONLY the continuous component; the atomic mass is handled elsewhere.
 */

package dr.inference.markovjumps;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

/**
 * @author Filippo Monti
 * @author Marc Suchard
 */

public class SericolaSeriesMarkovRewardFastModel extends AbstractModel {

    private static final boolean DEBUG = false;

    // Dependencies
    private final SubstitutionModel underlyingSubstitutionModel;

    // Dimensions
    private final int dim;
    private final int dim2;

    // -----------------------
    // Sorting buffers / maps
    // -----------------------
    private final double[] unsortedQ;     // dim2
    private final double[] sortedQ;       // dim2
    private final double[] sortedAlpha;   // dim

    // perm maps SORTED index -> ORIGINAL index
    private final int[] perm;
    // invPerm maps ORIGINAL index -> SORTED index
    private final int[] invPerm;

    // Mapping helpers for writing results directly to ORIGINAL order
    private final int[] outRowBaseBySorted; // perm[uS] * dim
    private final int[] outColBySorted;     // perm[vS]

    private boolean permDirty = true;
    private boolean qDirty = true;
    private boolean cachesDirty = true;

    // -----------------------
    // Numerics (mutable)
    // -----------------------
    private double lambda;                 // uniformization rate
    private final double[] P;              // dim2, P = I + sortedQ/lambda
    // invAlphaDiff[h] = 1/(sortedAlpha[h]-sortedAlpha[h-1]) for h=1..phi
    private final double[] invAlphaDiff;
    private final SericolaCumulantMatrices cumulantMatrices;
    private final SericolaRewardDensityPdf rewardDensityPdf;
    private final SericolaRewardDensityDerivative rewardDensityDerivative;
    private final ThreadLocal<SericolaRewardDensityWorkspace> workspaces;

    // Exponential for conditional probabilities (optional)
    private final EigenSystem eigenSystem;
    private EigenDecomposition eigenDecomposition;

    private final double[] out;

    private final boolean conditionalOnZ0;

    private final Parameter rewardRatesValues;   // distinct reward values
    private final Parameter rewardRatesValuesInternal;
    private final Parameter rewardRatesMapping;  // state i -> index into rewardRatesValues
    private final double[] alphaVals;
    private final int[] idx;

    public SericolaSeriesMarkovRewardFastModel(
            SubstitutionModel underlyingSubstitutionModel,
            Parameter rewardRatesValues,
            Parameter rewardRatesValuesInternal,
            Parameter rewardRatesMapping,
            int dim,
            double epsilon,
            boolean conditionalOnZ0
    ) {
        super("SericolaSeriesMarkovRewardFastModel");

        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("underlyingSubstitutionModel must be non-null");
        }
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesMapping == null) {
            throw new IllegalArgumentException("rewardRatesMapping must be non-null");
        }
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be > 0");
        }
        if (!(epsilon > 0.0 && epsilon < 1.0)) {
            throw new IllegalArgumentException("epsilon must be in (0,1)");
        }

        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesValuesInternal = rewardRatesValuesInternal;
        this.rewardRatesMapping = rewardRatesMapping;
        this.alphaVals = new double[dim];
        for (int i = 0; i < rewardRatesValues.getDimension(); i++) {
            final double v = rewardRatesValues.getParameterValue(i);
            if (v < 0 || v > 1) {
                throw new IllegalArgumentException("rewardRatesValues must be in [0,1], but value at index " + i + " is " + v);
            }
        }

        this.dim = dim;
        this.dim2 = dim * dim;

        this.unsortedQ = new double[dim2];
        this.sortedQ = new double[dim2];
        this.sortedAlpha = new double[dim];

        this.perm = new int[dim];
        this.invPerm = new int[dim];
        this.outRowBaseBySorted = new int[dim];
        this.outColBySorted = new int[dim];

        this.P = new double[dim2];
        this.invAlphaDiff = new double[dim];
        this.cumulantMatrices = new SericolaCumulantMatrices(dim);
        this.rewardDensityPdf =
                new SericolaRewardDensityPdf(dim, outRowBaseBySorted, outColBySorted);
        this.rewardDensityDerivative =
                new SericolaRewardDensityDerivative(dim, outRowBaseBySorted, outColBySorted);
        this.workspaces = ThreadLocal.withInitial(() -> new SericolaRewardDensityWorkspace(dim, epsilon));

        this.eigenSystem = new DefaultEigenSystem(dim);
        this.out = new double[dim2];
        this.conditionalOnZ0 = conditionalOnZ0; // TODO add this to the signature

        addModel(underlyingSubstitutionModel);
        addVariable(rewardRatesValues);
        addVariable(rewardRatesValuesInternal);
        addVariable(rewardRatesMapping);

        this.permDirty = true;
        this.qDirty = true;
        this.cachesDirty = true;
        this.idx = new int[dim];
        for (int i = 0; i < dim; i++) idx[i] = i;

    }

    // ============================================================
    // Public API: continuous density f^*(rewardProportion, t)
    // ============================================================

    public double[] computePdf(double rewardProportion, double time) {
        computePdfInto(rewardProportion, time, out);
        return out;
    }

    public void computePdfInto(double rewardProportion, double time, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        final SericolaRewardDensityWorkspace workspace = workspace();
        computePdfInto(
                workspace.scalarRewardProportions(rewardProportion),
                workspace.scalarTimes(time),
                false,
                workspace.scalarOutputs(out));
    }

    public void computePdfInto(double[] rewardProportions, double time, double[][] W) {
        computePdfInto(rewardProportions, workspace().scalarTimes(time), false, W);
    }

    public void computePdfInto(double[] rewardProportions, double[] times, double[][] W) {
        computePdfInto(rewardProportions, times, false, W);
    }

    public void computePdfInto(double[] rewardProportions, double[] times, boolean parsimonious, double[][] W) {
        validateInputs(rewardProportions, times, W);
        final int T = rewardProportions.length;
        if (T == 0) return;

        ensureNumericsUpToDate();
        validateAndZeroOutputs(W, T);

        final SericolaRewardDensityWorkspace workspace = workspace();
        final double maxT = workspace.preparePdf(
                rewardProportions,
                times,
                parsimonious,
                lambda,
                sortedAlpha,
                invAlphaDiff);

        // Ensure C is available up to required N (pdf uses n+1)
        ensureCForTime(maxT, /*extraN=*/1);
        final int N = cumulantMatrices.computedN() - 1;

        rewardDensityPdf.accumulateInto(
                W,
                T,
                N,
                lambda,
                invAlphaDiff,
                cumulantMatrices,
                workspace);

        if (conditionalOnZ0) {
            if (times.length == 1) {
                final double denom = -Math.expm1(-lambda * times[0]); // 1 - exp(-lambda t)
                if (denom <= 0.0) throw new IllegalStateException("Invalid denominator for conditional probability: " + denom);
                final double inv = 1.0 / denom;
                for (int t = 0; t < T; ++t) {
                    final double[] row = W[t];
                    for (int k = 0; k < dim2; ++k) row[k] *= inv;
                }
            } else {
                for (int t = 0; t < T; ++t) {
                    final double denom = -Math.expm1(-lambda * times[t]);
                    if (denom <= 0.0) throw new IllegalStateException("Invalid denominator for conditional probability at index " + t + ": " + denom);
                    final double inv = 1.0 / denom;
                    final double[] row = W[t];
                    for (int k = 0; k < dim2; ++k) row[k] *= inv;
                }
            }
        }
    }

    private void validateInputs(double[] rewardProportions, double[] times, double[][] W) {
        if (rewardProportions == null || times == null || W == null) {
            throw new IllegalArgumentException("rewardProportions/times/W must be non-null");
        }
        final int T = rewardProportions.length;
        final boolean singleTime = (times.length == 1);
        if (!singleTime && times.length != T) {
            throw new IllegalArgumentException("Either times.length==1 or times.length==rewardProportions.length");
        }
        if (W.length != T) {
            throw new IllegalArgumentException("W.length must equal rewardProportions.length");
        }
    }

    private void validateAndZeroOutputs(double[][] W, int T) {
        for (int t = 0; t < T; ++t) {
            final double[] row = W[t];
            if (row == null || row.length != dim2) {
                throw new IllegalArgumentException("W[" + t + "] must be non-null and length dim*dim=" + dim2);
            }
            Arrays.fill(row, 0.0);
        }
    }

    private void ensureCForTime(double time, int extraN) {
        final int requiredN = workspace().determineNumberOfSteps(lambda, time) + extraN;
        cumulantMatrices.ensureForTime(time, requiredN, P, sortedAlpha);
    }

    // ============================================================
    // Dirty-handling: permutation, sortedQ, numerics
    // ============================================================

    private void ensureNumericsUpToDate() {
        ensurePermutationUpToDate();
        ensureSortedQUpToDate();

        if (!cachesDirty) return;

        lambda = determineLambda(sortedQ);
        if (!(lambda > 0.0) || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalStateException("Invalid uniformization rate lambda=" + lambda);
        }
        fillP(sortedQ, lambda, P);

        Arrays.fill(invAlphaDiff, 0.0);
        for (int h = 1; h < dim; h++) {
            final double d = sortedAlpha[h] - sortedAlpha[h - 1];
            if (!(d > 0.0)) {
                throw new IllegalArgumentException("alphaRates must be strictly increasing after sorting; tie at h=" + h);
            }
            invAlphaDiff[h] = 1.0 / d;
        }

        eigenDecomposition = null;

        invalidateCComputedExtent();
        cachesDirty = false;
    }

    private void fillStateRewardRates(double[] stateRewardRates) {
        final double[] values = rewardRatesValues.getParameterValues();

        if (rewardRatesMapping.getDimension() != dim) {
            throw new IllegalArgumentException(
                    "rewardRatesMapping length mismatch: expected " + dim +
                            " got " + rewardRatesMapping.getDimension());
        }

        for (int i = 0; i < dim; i++) {
            final double raw = rewardRatesMapping.getParameterValue(i);
            final int map = (int) Math.round(raw);

            if (Math.abs(raw - map) > 1e-10) {
                throw new IllegalArgumentException(
                        "rewardRatesMapping must contain integer indices; found " + raw + " at state " + i);
            }
            if (map < 0 || map >= values.length) {
                throw new IllegalArgumentException(
                        "rewardRatesMapping[" + i + "]=" + map +
                                " is out of bounds for rewardRatesValues of length " + values.length);
            }

            final double a = values[map];
            if (a < 0.0 || a > 1.0 || Double.isNaN(a)) {
                throw new IllegalArgumentException(
                        "rewardRatesValues must lie in [0,1]; found " + a + " at value index " + map);
            }

            stateRewardRates[i] = a;
        }
    }

    private void ensurePermutationUpToDate() {
        if (!permDirty) return;

//        final double[] alphaVals = alphaRates.getParameterValues();
        fillStateRewardRates(alphaVals);

        if (alphaVals.length != dim) {
            throw new IllegalArgumentException("alphaRates length mismatch: expected " + dim + " got " + alphaVals.length);
        }

        for (int i = 0; i < dim; i++) {
            final double a = alphaVals[i];
            if (a < 0.0 || a > 1.0 || Double.isNaN(a)) {
                throw new IllegalArgumentException("alphaRates must lie in [0,1]; found " + a + " at index " + i);
            }
        }

//        Integer[] idx = new Integer[dim];
//        for (int i = 0; i < dim; i++) idx[i] = i;
//        Arrays.sort(idx, Comparator.comparingDouble(i -> alphaVals[i]));
        // idx already contains 0..dim-1 from the constructor
        sortIndicesByKey(idx, alphaVals);  // custom primitive sort

        for (int s = 0; s < dim; s++) {
            perm[s] = idx[s];
            sortedAlpha[s] = alphaVals[perm[s]];
        }
        for (int o = 0; o < dim; o++) {
            invPerm[o] = -1;
        }
        for (int s = 0; s < dim; s++) {
            invPerm[perm[s]] = s;
        }

        for (int uS = 0; uS < dim; uS++) {
            outRowBaseBySorted[uS] = perm[uS] * dim;
            outColBySorted[uS] = perm[uS];
        }

        permDirty = false;
        qDirty = true;
        cachesDirty = true;

        if (DEBUG) {
//            System.err.println("[DEBUG] alphaRates (original): " + Arrays.toString(alphaRates.getParameterValues()));
            System.err.println("[DEBUG] sortedAlpha: " + Arrays.toString(sortedAlpha));
            System.err.println("[DEBUG] perm (sorted -> original): " + Arrays.toString(perm));
            System.err.println("[DEBUG] invPerm (original -> sorted): " + Arrays.toString(invPerm));
        }
    }
    private static void sortIndicesByKey(int[] idx, double[] key) {
        quickSortByKey(idx, key, 0, idx.length - 1);
    }

    private static void quickSortByKey(int[] a, double[] key, int lo, int hi) {
        while (lo < hi) {
            int i = lo, j = hi;
            final double pivot = key[a[(lo + hi) >>> 1]];

            while (i <= j) {
                while (key[a[i]] < pivot) i++;
                while (key[a[j]] > pivot) j--;
                if (i <= j) {
                    final int tmp = a[i];
                    a[i] = a[j];
                    a[j] = tmp;
                    i++;
                    j--;
                }
            }

            // Recurse on smaller partition first (limits stack depth)
            if (j - lo < hi - i) {
                if (lo < j) quickSortByKey(a, key, lo, j);
                lo = i;
            } else {
                if (i < hi) quickSortByKey(a, key, i, hi);
                hi = j;
            }
        }
    }

    private void ensureSortedQUpToDate() {
        if (!qDirty) return;

        underlyingSubstitutionModel.getInfinitesimalMatrix(unsortedQ);
        if (unsortedQ.length != dim2) {
            throw new IllegalStateException("Internal unsortedQ length mismatch");
        }

        for (int iS = 0; iS < dim; ++iS) {
            final int iO = perm[iS];
            final int rowS = iS * dim;
            final int rowO = iO * dim;
            for (int jS = 0; jS < dim; ++jS) {
                final int jO = perm[jS];
                sortedQ[rowS + jS] = unsortedQ[rowO + jO];
            }
        }

        qDirty = false;
        cachesDirty = true;
    }

    // ============================================================
    // Helpers
    // ============================================================

    private SericolaRewardDensityWorkspace workspace() {
        return workspaces.get();
    }

    private double determineLambda(double[] QrowMajor) {
        double minDiag = QrowMajor[0];
        for (int i = 1; i < dim; ++i) {
            final double d = QrowMajor[i * dim + i];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
    }

    private void invalidateCComputedExtent() {
        cumulantMatrices.invalidateComputedExtent();
    }

    private void fillP(double[] Qsorted, double lambda, double[] Pout) {
        final double inv = 1.0 / lambda;
        for (int i = 0; i < dim; ++i) {
            final int ioff = i * dim;
            for (int j = 0; j < dim; ++j) {
                Pout[ioff + j] = (i == j ? 1.0 : 0.0) + Qsorted[ioff + j] * inv;
            }
        }
    }

    public double getLambda() {
        ensureNumericsUpToDate();
        return lambda;
    }

    public double getUniformizationRate() {
        return getLambda();
    }

    private double[][] squareMatrix(final double[] matRowMajor) {
        final double[][] rtn = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(matRowMajor, i * dim, rtn[i], 0, dim);
        }
        return rtn;
    }

    // ============================================================
    // Optional: conditional probabilities exp(Q t)
    // ============================================================

    private EigenDecomposition getEigenDecomposition() {
        ensureNumericsUpToDate();
        if (eigenDecomposition == null) {
            eigenDecomposition = eigenSystem.decomposeMatrix(squareMatrix(sortedQ));
        }
        return eigenDecomposition;
    }

    public double[] computeConditionalProbabilities(double distance) {
        ensureNumericsUpToDate();
        final double[] matrix = new double[dim2];
        eigenSystem.computeExponential(getEigenDecomposition(), distance, matrix);
        return matrix;
    }

    public double computeConditionalProbability(double distance, int i, int j) {
        ensureNumericsUpToDate();
        return eigenSystem.computeExponential(getEigenDecomposition(), distance, i, j);
    }

    // ============================================================
    // BEAST Model callbacks
    // ============================================================

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == underlyingSubstitutionModel) {
            qDirty = true;
            cachesDirty = true;
            fireModelChanged();
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rewardRatesMapping || variable == rewardRatesValues || variable == rewardRatesValuesInternal) {
            permDirty = true;
            qDirty = true;
            cachesDirty = true;
            fireModelChanged();
        }
    }

    @Override
    protected void storeState() {
        // nothing needed (we recompute on restore)
    }

    @Override
    protected void restoreState() {
        permDirty = true;
        qDirty = true;
        cachesDirty = true;
        invalidateCComputedExtent();
        eigenDecomposition = null;
    }

    @Override
    protected void acceptState() {
        // nothing
    }

    @Override
    public String toString() {
        ensureNumericsUpToDate();
        StringBuilder sb = new StringBuilder();
        sb.append("sortedQ: ").append(new Vector(sortedQ)).append("\n");
        sb.append("sortedAlpha: ").append(new Vector(sortedAlpha)).append("\n");
        sb.append("lambda: ").append(lambda).append("\n");
        sb.append("allocatedN: ").append(cumulantMatrices.allocatedN()).append("\n");
        final double maxTime = cumulantMatrices.maxTime();
        sb.append("maxTime: ").append(maxTime).append("\n");
        sb.append("cprob at maxTime: ").append(new Vector(computeConditionalProbabilities(maxTime))).append("\n");
        return sb.toString();
    }

    public void computePdfDerivativeWrtYInto(double rewardProportion, double time, double[] differential) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    // ============================================================
    // NO-ALLOCATION derivative of reward-density w.r.t. rewardProportion
    // ============================================================
    //
    // Computes: d/drewardProportion f^*_{ij}(rewardProportion | t)
    //
    // One-sided interior derivative code exists in SericolaRewardDensityDerivative,
    // but exact support-boundary derivative calls are rejected until that convention
    // is made explicit for the full reward-aware gradient.

    public void computePdfDerivativeWrtRewardProportionInto(double rewardProportion, double branchLength, double[] out, boolean shiftback) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        if (branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be > 0");
        }
        Arrays.fill(out, 0.0);
        ensureNumericsUpToDate();
        final SericolaRewardDensityWorkspace workspace = workspace();
        final int h = workspace.prepareDerivative(rewardProportion, sortedAlpha, invAlphaDiff);
        if (workspace.isZero(0) || workspace.isOne(0)) {
            throw new UnsupportedOperationException(
                    "Boundary Values for rewardProportion touched; one-sided interior derivative is implemented but disabled " +
                            "until the boundary convention is explicit. rewardProportion=" + rewardProportion +
                            ", interval=[" + sortedAlpha[h - 1] + "," + sortedAlpha[h] + "], h=" + h +
                            ", xh=" + workspace.xh(0));
        }

        ensureCForTime(branchLength, /*extraN=*/1);
        final int N = cumulantMatrices.computedN() - 1;

        rewardDensityDerivative.computeWrtRewardProportionInto(
                h,
                N,
                branchLength,
                lambda,
                invAlphaDiff[h],
                workspace.xh(0),
                workspace.isZero(0),
                workspace.isOne(0),
                cumulantMatrices,
                workspace.increments(),
                out);
    }

}
