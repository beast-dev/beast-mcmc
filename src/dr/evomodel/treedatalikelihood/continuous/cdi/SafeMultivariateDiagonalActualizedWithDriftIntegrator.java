package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class SafeMultivariateDiagonalActualizedWithDriftIntegrator extends SafeMultivariateWithDriftIntegrator {

    private static boolean DEBUG = false;

    public SafeMultivariateDiagonalActualizedWithDriftIntegrator(PrecisionType precisionType,
                                                                 int numTraits, int dimTrait, int dimProcess,
                                                                 int bufferCount, int diffusionCount) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateDiagonalActualizedWithDriftIntegrator");
    }

    @Override
    public void getBranch1mActualization(int bufferIndex, double[] diagonal1mActualization) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (diagonal1mActualization != null);
        assert (diagonal1mActualization.length >= dimTrait);

        System.arraycopy(diagonal1mActualizations, bufferIndex * dimTrait,
                diagonal1mActualization, 0, dimTrait);
    }

    @Override
    public void getBranchActualization(int bufferIndex, double[] diagonalActualization) {
        getBranch1mActualization(bufferIndex, diagonalActualization);
        oneMinus(diagonalActualization);
    }

    @Override
    public void getBranchExpectation(double[] actualization, double[] parentValue, double[] displacement,
                                     double[] expectation) {

        assert (expectation != null);
        assert (expectation.length >= dimTrait);

        assert (actualization != null);
        assert (actualization.length >= dimTrait);

        assert (parentValue != null);
        assert (parentValue.length >= dimTrait);

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        for (int i = 0; i < dimTrait; ++i) {
            expectation[i] = actualization[i] * parentValue[i] + displacement[i];
        }
    }

    private static final boolean TIMING = false;

    private void allocateStorage() {

        diagonal1mActualizations = new double[dimTrait * bufferCount];
        stationaryVariances = new double[dimProcess * dimProcess * diffusionCount];

        vectorDiagQdi = new double[dimTrait];
        vectorDiagQdj = new double[dimTrait];
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Setting variances, displacement and actualization vectors
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex,
                                               final double[] diagonalAlpha, final double[] rotation) {

        assert (stationaryVariances != null);

        assert dimProcess == diagonalAlpha.length;

        final int matrixSize = dimProcess * dimProcess;
        final int offset = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(diagonalAlpha, scales);

        setStationaryVariance(offset, scales, matrixSize, rotation);
    }

    void setStationaryVariance(int offset, double[] scales, int matrixSize, double[] rotation) {
        scaleInv(inverseDiffusions, offset, scales, stationaryVariances, offset, matrixSize);
    }

    static void scaleInv(final double[] source,
                         final int sourceOffset,
                         final double[] scales,
                         final double[] destination,
                         final int destinationOffset,
                         final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] = source[sourceOffset + i] / scales[i];
        }
    }

    private static void scalingMatrix(double[] eigAlpha, double[] scales) {
        int nEig = eigAlpha.length;
        for (int i = 0; i < nEig; ++i) {
            for (int j = 0; j < nEig; ++j) {
                scales[i * nEig + j] = eigAlpha[i] + eigAlpha[j];
            }
        }
    }

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                         final double[] edgeLengths, final double[] optimalRates,
                                                         final double[] diagonalStrengthOfSelectionMatrix,
                                                         final double[] rotation,
                                                         int updateCount) {

        assert (diffusions != null);
        assert (probabilityIndices.length >= updateCount);
        assert (edgeLengths.length >= updateCount);

        super.updateOrnsteinUhlenbeckDiffusionMatrices(precisionIndex, probabilityIndices, edgeLengths, optimalRates,
                diagonalStrengthOfSelectionMatrix, rotation, updateCount);

        if (DEBUG) {
            System.err.println("Matrices (safe with actualized drift):");
        }

        final int matrixTraitSize = dimTrait * dimTrait;
        final int matrixProcessSize = dimProcess * dimProcess;
        final int unscaledOffset = matrixProcessSize * precisionIndex;


        if (TIMING) {
            startTime("actualization");
        }

        for (int up = 0; up < updateCount; ++up) {

            final double edgeLength = edgeLengths[up];

            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];
            final int scaledOffset = dimTrait * scaledOffsetDiagonal;

            computeOUActualization(diagonalStrengthOfSelectionMatrix, rotation, edgeLength,
                    scaledOffsetDiagonal, scaledOffset);
        }

        if (TIMING) {
            endTime("actualization");
        }


        if (TIMING) {
            startTime("diffusion");
        }

        for (int up = 0; up < updateCount; ++up) {

            final double edgeLength = edgeLengths[up];

            final int scaledOffset = matrixTraitSize * probabilityIndices[up];
            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];

            computeOUVarianceBranch(unscaledOffset, scaledOffset, scaledOffsetDiagonal, edgeLength);

            invertVectorSymmPosDef(variances, precisions, scaledOffset, dimProcess);
        }

        if (TIMING) {
            endTime("diffusion");
        }

        assert (optimalRates != null);
        assert (displacements != null);
        assert (optimalRates.length >= updateCount * dimProcess);

        if (TIMING) {
            startTime("drift1");
        }

        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {

            final int pio = dimTrait * probabilityIndices[up];
            final int scaledOffset = matrixTraitSize * probabilityIndices[up];

            computeOUActualizedDisplacement(optimalRates, offset, scaledOffset, pio);
            offset += dimProcess;
        }

        if (TIMING) {
            endTime("drift1");
        }
    }

    void computeOUActualization(final double[] diagonalStrengthOfSelectionMatrix,
                                final double[] rotation,
                                final double edgeLength,
                                final int scaledOffsetDiagonal,
                                final int scaledOffset) {
        computeOUDiagonal1mActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait,
                diagonal1mActualizations, scaledOffsetDiagonal);
    }

    static void computeOUDiagonal1mActualization(final double[] source,
                                               final double edgeLength,
                                               final int dim,
                                               final double[] destination,
                                               final int destinationOffset) {
        for (int i = 0; i < dim; ++i) {
            destination[destinationOffset + i] = - Math.expm1(-source[i] * edgeLength);
        }
    }

    void computeOUVarianceBranch(final int unscaledOffset,
                               final int scaledOffset,
                               final int scaledOffsetDiagonal,
                               final double edgeLength) {
        scalingActualizationMatrix(diagonal1mActualizations, scaledOffsetDiagonal,
                stationaryVariances, unscaledOffset,
                variances, scaledOffset,
                dimTrait,
                edgeLength,
                inverseDiffusions,
                unscaledOffset);
    }

    private static void scalingActualizationMatrix(final double[] diagonal1mActualizations,
                                                   final int offsetActualization,
                                                   final double[] stationaryVariances,
                                                   final int offsetStationaryVariances,
                                                   final double[] destination,
                                                   final int destinationOffset,
                                                   final int dim,
                                                   final double edgeLength,
                                                   final double[] variance,
                                                   final int varianceOffset) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double var = stationaryVariances[offsetStationaryVariances + i * dim + j];
                if (Double.isInfinite(var) || (diagonal1mActualizations[offsetActualization + i] + diagonal1mActualizations[offsetActualization + j] == 0.0)) {
                    destination[destinationOffset + i * dim + j] = edgeLength * variance[varianceOffset + i * dim + j];
                } else {
                    destination[destinationOffset + i * dim + j] = var * (- diagonal1mActualizations[offsetActualization + i] * diagonal1mActualizations[offsetActualization + j] + diagonal1mActualizations[offsetActualization + i] + diagonal1mActualizations[offsetActualization + j]);
                }
            }
        }
    }

    private static void invertVectorSymmPosDef(final double[] source,
                                               final double[] destination,
                                               final int offset,
                                               final int dim) {
        DenseMatrix64F sourceMatrix = wrap(source, offset, dim, dim);
        DenseMatrix64F destinationMatrix = new DenseMatrix64F(dim, dim);

//        CommonOps.invert(sourceMatrix, destinationMatrix);
        symmPosDefInvert(sourceMatrix, destinationMatrix);

        unwrap(destinationMatrix, destination, offset);
    }

    void computeOUActualizedDisplacement(final double[] optimalRates,
                                         final int offset,
                                         final int actualizationOffset,
                                         final int pio) {

        for (int j = 0; j < dimTrait; ++j) {
            displacements[pio + j] = optimalRates[offset + j] * diagonal1mActualizations[pio + j];
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Tree-traversal functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    void actualizePrecision(DenseMatrix64F Pjp, DenseMatrix64F QjPjp, int jbo, int jmo, int jdo) {
        final double[] diagQdj = vectorDiagQdj;
        System.arraycopy(diagonal1mActualizations, jdo, diagQdj, 0, dimTrait);
        oneMinus(diagQdj);
        MissingOps.diagMult(diagQdj, Pjp, QjPjp);
        MissingOps.diagMult(QjPjp, diagQdj, Pjp);
    }

    @Override
    void actualizeVariance(DenseMatrix64F Vip, int ibo, int imo, int ido) {
        final double[] diagQdi = vectorDiagQdi;
        System.arraycopy(diagonal1mActualizations, ido, diagQdi, 0, dimTrait);
        oneMinus(diagQdi);
        diagonalDoubleProduct(Vip, diagQdi, Vip);
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] = (1.0 - diagonal1mActualizations[ido + g]) * preOrderPartials[ibo + g] + displacements[ido + g];
        }
    }

    public double[] getStationaryVariance(int precisionIndex) {

        assert (stationaryVariances != null);

        return getMatrixProcess(precisionIndex, stationaryVariances);
    }

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {

        final double[] diagQdi = vectorDiagQdi;
        System.arraycopy(diagonal1mActualizations, ido, diagQdi, 0, dimTrait);
        oneMinus(diagQdi);
        final double[] diagQdj = vectorDiagQdj;
        System.arraycopy(diagonal1mActualizations, jdo, diagQdj, 0, dimTrait);
        oneMinus(diagQdj);

        final DenseMatrix64F QdiPipQdi = matrix0;
        final DenseMatrix64F QdjPjpQdj = matrix1;
        diagonalDoubleProduct(Pip, diagQdi, QdiPipQdi);
        diagonalDoubleProduct(Pjp, diagQdj, QdjPjpQdj);
        CommonOps.add(QdiPipQdi, QdjPjpQdj, Pk);

        if (DEBUG) {
            System.err.println("Qdi: " + Arrays.toString(diagQdi));
            System.err.println("\tQdiPipQdi: " + QdiPipQdi);
            System.err.println("\tQdj: " + Arrays.toString(diagQdj));
            System.err.println("\tQdjPjpQdj: " + QdjPjpQdj);
        }
    }

    @Override
    void computeWeightedSum(final double[] ipartial,
                            final double[] jpartial,
                            final int dimTrait,
                            final double[] out) {
        weightedSumActualized(ipartial, 0, matrixPip, vectorDiagQdi, 0,
                jpartial, 0, matrixPjp, vectorDiagQdj, 0,
                dimTrait, out);
    }

    private static void diagonalDoubleProduct(DenseMatrix64F source, double[] diagonalScales,
                                              DenseMatrix64F destination) {
        MissingOps.diagMult(source, diagonalScales, destination);
        MissingOps.diagMult(diagonalScales, destination);
    }

    static void oneMinus(double[] x) {
        for (int i = 0; i < x.length; i++) {
            x[i] = 1.0 - x[i];
        }
    }

    private double[] diagonal1mActualizations;
    double[] stationaryVariances;
    private double[] vectorDiagQdi;
    private double[] vectorDiagQdj;
}
