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
                                                                 int numTraits, int dimTrait, int dimProcess,
                                                                 int bufferCount, int diffusionCount) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

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
        computeOUDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait,
                diagonalActualizations, scaledOffsetDiagonal);
    }

    static void computeOUDiagonalActualization(final double[] source,
                                               final double edgeLength,
                                               final int dim,
                                               final double[] destination,
                                               final int destinationOffset) {
        for (int i = 0; i < dim; ++i) {
            destination[destinationOffset + i] = Math.exp(-source[i] * edgeLength);
        }
    }

    void computeOUVarianceBranch(final int unscaledOffset,
                               final int scaledOffset,
                               final int scaledOffsetDiagonal,
                               final double edgeLength) {
        scalingActualizationMatrix(diagonalActualizations, scaledOffsetDiagonal,
                stationaryVariances, unscaledOffset,
                variances, scaledOffset,
                dimTrait,
                edgeLength,
                inverseDiffusions,
                unscaledOffset);
    }

    private static void scalingActualizationMatrix(final double[] diagonalActualizations,
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
                if (Double.isInfinite(var) || (1 - diagonalActualizations[offsetActualization + i] * diagonalActualizations[offsetActualization + j]) == 0.0) {
                    destination[destinationOffset + i * dim + j] = edgeLength * variance[varianceOffset + i * dim + j];
                } else {
                    destination[destinationOffset + i * dim + j] = var * (1 - diagonalActualizations[offsetActualization + i] * diagonalActualizations[offsetActualization + j]);
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
