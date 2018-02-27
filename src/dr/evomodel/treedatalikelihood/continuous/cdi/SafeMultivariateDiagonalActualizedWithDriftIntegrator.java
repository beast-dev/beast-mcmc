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

    public SafeMultivariateDiagonalActualizedWithDriftIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
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

    private static final boolean TIMING = false;

    private void allocateStorage() {

        diagonalActualizations = new double[dimTrait * bufferCount];
        stationaryVariances = new double[dimTrait * dimTrait * diffusionCount];

        vectorDiagQdi = new double[dimTrait];
        vectorDiagQdj = new double[dimTrait];
    }


    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] diagonalAlpha, final double[] rotation) {

        assert (stationaryVariances != null);

        final int matrixSize = dimTrait * dimTrait;
        final int offset = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(diagonalAlpha, scales);

        if (rotation.length == 0) {
            scaleInv(inverseDiffusions, offset, scales, stationaryVariances, offset, matrixSize);
        } else {
            assert (rotation.length == matrixSize);
            DenseMatrix64F rotMat = wrap(rotation, 0, dimTrait, dimTrait);
            DenseMatrix64F variance = wrap(inverseDiffusions, offset, dimTrait, dimTrait);
            transformMatrix(variance, rotMat);
            double[] transVar = new double[matrixSize];
            unwrap(variance, transVar, 0);
            scaleInv(transVar, 0, scales, stationaryVariances, offset, matrixSize);
        }
    }

//    private void transformMatrix(DenseMatrix64F matrix, DenseMatrix64F rotation) {
//        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
//        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
//        CommonOps.mult(rotation, matrix, tmp);
//        CommonOps.multTransB(tmp, rotation, matrix);
//    }

    private void transformMatrix(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransA(rotation, matrix, tmp); // TODO: this is a specialized version for symmetric A (see above general version)
        CommonOps.mult(tmp, rotation, matrix);
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

    private static void scaleInv(final double[] source,
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

        if (DEBUG) {
            System.err.println("Matrices (safe with drift):");
        }

        final int matrixSize = dimTrait * dimTrait;
        final int unscaledOffset = matrixSize * precisionIndex;


        if (TIMING) {
            startTime("actualization");
        }

        for (int up = 0; up < updateCount; ++up) {

            if (DEBUG) {
                System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
            }

            final double edgeLength = edgeLengths[up];
            branchLengths[dimMatrix * probabilityIndices[up]] = edgeLength;  // TODO Remove dimMatrix

            final int scaledOffset = dimTrait * probabilityIndices[up];

            computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait, diagonalActualizations, scaledOffset);
        }

        if (TIMING) {
            endTime("actualization");
        }


        if (TIMING) {
            startTime("diffusion");
        }

        for (int up = 0; up < updateCount; ++up) {

            if (DEBUG) {
                System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
            }

            final int scaledOffset = matrixSize * probabilityIndices[up];
            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];

            double[] scales = new double[matrixSize];
            scalingActualizationMatrix(diagonalActualizations, scaledOffsetDiagonal, dimTrait, scales);

            scale(stationaryVariances, unscaledOffset, scales, variances, scaledOffset, matrixSize);
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


        for (int up = 0; up < updateCount; ++up) {
            final int pio = dimTrait * probabilityIndices[up];
            for (int j = 0; j < dimTrait; ++j) {
                displacements[pio + j] = optimalRates[up + j] * (1 - diagonalActualizations[pio + j]);
            }
        }

        if (TIMING) {
            endTime("drift1");
        }

        // NOTE TO PB: very complex function, why multiple for (up = 0; up < updateCount; ++up) ?
        // PB: I don't know if that's better, but I needed the loop to re-order the optimal rates according to the branches.
        // But maybe there is a better way to do that ?

        precisionOffset = dimTrait * dimTrait * precisionIndex;
        precisionLogDet = determinants[precisionIndex];
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


//    private static void scalingActualizationDisplacement(final double[] vec) {
//        for (int i = 0; i < vec.length; ++i) {
//            vec[i] = 1 - vec[i];
//        }
//    }

//    private static void scalingActualizationDisplacement(final double[] source,
//                                                         final int offset,
//                                                         final int dim,
//                                                         final double[] destination) {
//        for (int i = 0; i < dim; ++i) {
//            destination[i] = 1 - source[offset + i];
//        }
//    }

    @Override
    public void updatePreOrderPartial(
            final int kBuffer, // parent
            final int iBuffer, // node
            final int iMatrix,
            final int jBuffer, // sibling
            final int jMatrix) {

        throw new RuntimeException("Not yet implemented");
    }

    private static void diagonalDoubleProduct(DenseMatrix64F source, double[] diagonalScales, DenseMatrix64F destination) {
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
                scales[i * dim + j] = diagonalScales[i] * diagonalScales[j];
            }
        }
    }

    private double[] diagonalActualizations;
    double[] stationaryVariances;
    private double[] vectorDiagQdi;
    private double[] vectorDiagQdj;

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
}
