package dr.evomodel.treedatalikelihood.continuous.cdi;

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
                                                                 int numTraits, int dimTrait, int bufferCount,
                                                                 int diffusionCount) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateDiagonalActualizedWithDriftIntegrator");
    }

    // NOTE TO PB: need to merge all SafeMultivariate* together and then delegate specialized work ... avoid massive code duplication
    // PB: Merged `updatePartial` functions by breaking it into smaller functions.


    @Override
    public void getBranchActualization(int bufferIndex, double[] diagonalActualization) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (diagonalActualization != null);
        assert (diagonalActualization.length >= dimTrait);

        System.arraycopy(diagonalActualizations, bufferIndex * dimTrait,
                diagonalActualization, 0, dimTrait);
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

        diagonalActualizations = new double[dimTrait * bufferCount];
        stationaryVariances = new double[dimTrait * dimTrait * diffusionCount];

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

        final int matrixSize = dimTrait * dimTrait;
        final int offset = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(diagonalAlpha, scales);

        setStationaryVariance(offset, scales, matrixSize, rotation);
    }

    void setStationaryVariance(int offset, double[] scales, int matrixSize, double[] rotation) {
        scaleInv(inverseDiffusions, offset, scales, stationaryVariances, offset, matrixSize);
    }

    private static void scale(final double[] source,
                              final int sourceOffset,
                              final double[] scales,
                              final double[] destination,
                              final int destinationOffset,
                              final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] = scales[i] * source[sourceOffset + i];
        }
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

        final int matrixSize = dimTrait * dimTrait;
        final int unscaledOffset = matrixSize * precisionIndex;


        if (TIMING) {
            startTime("actualization");
        }

        for (int up = 0; up < updateCount; ++up) {

            final double edgeLength = edgeLengths[up];

            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];
            final int scaledOffset = dimTrait * scaledOffsetDiagonal;

            computeActualization(diagonalStrengthOfSelectionMatrix, rotation, edgeLength,
                    scaledOffsetDiagonal, scaledOffset);
        }

        if (TIMING) {
            endTime("actualization");
        }


        if (TIMING) {
            startTime("diffusion");
        }

        for (int up = 0; up < updateCount; ++up) {

            final int scaledOffset = matrixSize * probabilityIndices[up];
            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];

            computeVarianceBranch(unscaledOffset, scaledOffset, scaledOffsetDiagonal);

            invertVector(variances, precisions, scaledOffset, dimTrait);
        }

        if (TIMING) {
            endTime("diffusion");
        }

        assert (optimalRates != null);
        assert (displacements != null);
        assert (optimalRates.length >= updateCount * dimTrait);

        if (TIMING) {
            startTime("drift1");
        }

//        int offset = 0;
//        for (int up = 0; up < updateCount; ++up) {
//
//            final int pio = dimTrait * probabilityIndices[up];
//
//            double[] scales = new double[dimTrait];
//            scalingActualizationDisplacement(diagonalActualizations, pio, dimTrait, scales);
//
//            scale(optimalRates, offset, scales, displacements, pio, dimTrait);
//            offset += dimTrait;
//        }

        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {

            final int pio = dimTrait * probabilityIndices[up];
            final int scaledOffset = matrixSize * probabilityIndices[up];

            computeActualizedDisplacement(optimalRates, offset, scaledOffset, pio);
            offset += dimTrait;
        }

        if (TIMING) {
            endTime("drift1");
        }

        // NOTE TO PB: very complex function, why multiple for (up = 0; up < updateCount; ++up) ?
        // PB: I don't know if that's better, but I needed the loop to re-order the optimal rates according to the branches.
        // Also now re-factored to avoid code duplication in SafeMultivariateActualized*
    }

    void computeActualization(final double[] diagonalStrengthOfSelectionMatrix,
                              final double[] rotation,
                              final double edgeLength,
                              final int scaledOffsetDiagonal,
                              final int scaledOffset) {
        computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait,
                diagonalActualizations, scaledOffsetDiagonal);
    }

    static void computeDiagonalActualization(final double[] source,
                                             final double edgeLength,
                                             final int dim,
                                             final double[] destination,
                                             final int destinationOffset) {
        for (int i = 0; i < dim; ++i) {
            destination[destinationOffset + i] = Math.exp(-source[i] * edgeLength);
        }
    }

    void computeVarianceBranch(final int sourceOffset,
                               final int destinationOffset,
                               final int destinationOffsetDiagonal) {
        double[] scales = new double[dimTrait * dimTrait];
        scalingActualizationMatrix(diagonalActualizations, destinationOffsetDiagonal, dimTrait, scales);
        scale(stationaryVariances, sourceOffset, scales, variances, destinationOffset, dimTrait * dimTrait);
    }

    private static void scalingActualizationMatrix(final double[] source,
                                                   final int offset,
                                                   final int dim,
                                                   final double[] destination) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                destination[i * dim + j] = 1 - source[offset + i] * source[offset + j];
            }
        }
    }

    private static void invertVector(final double[] source,
                                     final double[] destination,
                                     final int offset,
                                     final int dim) {
        DenseMatrix64F sourceMatrix = wrap(source, offset, dim, dim);
        DenseMatrix64F destinationMatrix = new DenseMatrix64F(dim, dim);

        CommonOps.invert(sourceMatrix, destinationMatrix);

        unwrap(destinationMatrix, destination, offset);
    }

    void computeActualizedDisplacement(final double[] optimalRates,
                                       final int offset,
                                       final int actualizationOffset,
                                       final int pio) {

        for (int j = 0; j < dimTrait; ++j) {
            displacements[pio + j] = optimalRates[offset + j] * (1 - diagonalActualizations[pio + j]);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Tree-traversal functions
    ///////////////////////////////////////////////////////////////////////////

//    @Override
//    public void updatePreOrderPartial(
//            final int kBuffer, // parent
//            final int iBuffer, // node
//            final int iMatrix,
//            final int jBuffer, // sibling
//            final int jMatrix) {
//
//        throw new RuntimeException("Not yet implemented");
//    }

    @Override
    void actualizePrecision(DenseMatrix64F Pjp, DenseMatrix64F QjPjp, int jbo, int jmo, int jdo) {
        final double[] diagQdj = vectorDiagQdj;
        System.arraycopy(diagonalActualizations, jdo, diagQdj, 0, dimTrait);
        diagonalProduct(Pjp, diagQdj, QjPjp);
        diagonalDoubleProduct(Pjp, diagQdj, Pjp);
    }

    @Override
    void actualizeVariance(DenseMatrix64F Vip, int ibo, int imo, int ido) {
        final double[] diagQdi = vectorDiagQdi;
        System.arraycopy(diagonalActualizations, ido, diagQdi, 0, dimTrait);
        diagonalDoubleProduct(Vip, diagQdi, Vip);
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] = diagonalActualizations[ido + g] * preOrderPartials[ibo + g] + displacements[ido + g];
        }
    }

    public double[] getStationaryVariance(int precisionIndex) {

        assert (stationaryVariances != null);

        final int offset = dimTrait * dimTrait * precisionIndex;

        double[] buffer = new double[dimTrait * dimTrait];

        System.arraycopy(stationaryVariances, offset, buffer, 0, dimTrait * dimTrait);

        return buffer;
    }

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {

        final double[] diagQdi = vectorDiagQdi;
        System.arraycopy(diagonalActualizations, ido, diagQdi, 0, dimTrait);
        final double[] diagQdj = vectorDiagQdj;
        System.arraycopy(diagonalActualizations, jdo, diagQdj, 0, dimTrait);

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
        double[] scales = new double[diagonalScales.length * diagonalScales.length];
        diagonalToMatrixDouble(diagonalScales, scales);

        for (int i = 0; i < destination.data.length; ++i) {
            destination.data[i] = scales[i] * source.data[i];
        }
    }

    private static void diagonalToMatrixDouble(double[] diagonalScales, double[] scales) {
        int dim = diagonalScales.length;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                scales[i * dim + j] = diagonalScales[i] * diagonalScales[j];
            }
        }
    }

    private static void diagonalProduct(DenseMatrix64F source, double[] diagonalScales,
                                        DenseMatrix64F destination) {
        double[] scales = new double[diagonalScales.length * diagonalScales.length];
        diagonalToMatrix(diagonalScales, scales);

        for (int i = 0; i < destination.data.length; ++i) {
            destination.data[i] = scales[i] * source.data[i];
        }
    }

    private static void diagonalToMatrix(double[] diagonalScales, double[] scales) {
        int dim = diagonalScales.length;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                scales[i * dim + j] = diagonalScales[i];
            }
        }
    }

    private double[] diagonalActualizations;
    double[] stationaryVariances;
    private double[] vectorDiagQdi;
    private double[] vectorDiagQdj;
}
