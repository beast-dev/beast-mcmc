/*
 * SafeMultivariateDiagonalActualizedWithDriftIntegrator.java
 *
 * Diagonal specialisation of SafeMultivariateActualizedWithDriftIntegrator:
 * Q(t) is diagonal in trait space.
 *
 * Copyright © 2002-2024 the BEAST Development Team
 */

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
public class SafeMultivariateDiagonalActualizedWithDriftIntegrator
        extends SafeMultivariateActualizedWithDriftIntegrator {

    private static final boolean DEBUG = false;

    /** (1 - exp(-λ_i t)) per trait and branch. */
    private double[] diagonal1mActualizations; // dimTrait x bufferCount

    /** Temporary vectors for tree traversal. */
    private double[] vectorDiagQdi;
    private double[] vectorDiagQdj;

    public SafeMultivariateDiagonalActualizedWithDriftIntegrator(
            PrecisionType precisionType,
            int numTraits,
            int dimTrait,
            int dimProcess,
            int bufferCount,
            int diffusionCount) {

        super(precisionType, numTraits, dimTrait, dimProcess,
                bufferCount, diffusionCount,
                true /* symmetric basis */);

        allocateDiagonalStorage();
        // System.err.println("Trying SafeMultivariateDiagonalActualizedWithDriftIntegrator");
    }

    private void allocateDiagonalStorage() {
        diagonal1mActualizations = new double[dimTrait * bufferCount];
        vectorDiagQdi            = new double[dimTrait];
        vectorDiagQdj            = new double[dimTrait];
    }

    /* **********************************************************************
     * Branch accessors
     * ******************************************************************* */

    @Override
    public void getBranch1mActualization(int bufferIndex,
                                         double[] diagonal1mActualization) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert diagonal1mActualization != null;
        assert diagonal1mActualization.length >= dimTrait;

        System.arraycopy(diagonal1mActualizations,
                bufferIndex * dimTrait,
                diagonal1mActualization, 0,
                dimTrait);
    }

    @Override
    public void getBranchActualization(int bufferIndex,
                                       double[] diagonalActualization) {
        getBranch1mActualization(bufferIndex, diagonalActualization);
        oneMinus(diagonalActualization);
    }

    @Override
    public void getBranchExpectation(double[] actualization,
                                     double[] parentValue,
                                     double[] displacement,
                                     double[] expectation) {

        assert expectation   != null && expectation.length   >= dimTrait;
        assert actualization != null && actualization.length >= dimTrait;
        assert parentValue   != null && parentValue.length   >= dimTrait;
        assert displacement  != null && displacement.length  >= dimTrait;

        for (int i = 0; i < dimTrait; ++i) {
            expectation[i] = actualization[i] * parentValue[i] + displacement[i];
        }
    }

    /* **********************************************************************
     * Stationary variance Σ – diagonal variant
     * ******************************************************************* */

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex,
                                               final double[] diagonalAlpha,
                                               final double[] rotation) {

        assert stationaryVariances != null;
        assert dimProcess == diagonalAlpha.length;

        final int matrixSize = dimProcess * dimProcess;
        final int offset     = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(diagonalAlpha, scales);

        // In the diagonal case, the precision is already in the correct basis:
        setStationaryVariance(offset, scales, matrixSize, rotation);
    }

    @Override
    void setStationaryVariance(int offset,
                               double[] scales,
                               int matrixSize,
                               double[] rotation) {
        // Σ = invDiffusion / (λ_i + λ_j) directly, no basis change.
        scaleInv(inverseDiffusions, offset,
                scales,
                stationaryVariances, offset,
                matrixSize);
    }

    /* **********************************************************************
     * OU hooks – diagonal specialisation
     * ******************************************************************* */

    @Override
    void computeOUActualization(final double[] diagonalStrengthOfSelectionMatrix,
                                final double[] rotation,
                                final double edgeLength,
                                final int scaledOffsetDiagonal,
                                final int scaledOffset) {

        computeOUDiagonal1mActualization(diagonalStrengthOfSelectionMatrix,
                edgeLength,
                dimTrait,
                diagonal1mActualizations,
                scaledOffsetDiagonal);
    }

    @Override
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

    @Override
    void computeOUActualizedDisplacement(final double[] optimalRates,
                                         final int offset,
                                         final int actualizationOffset,
                                         final int pio) {

        for (int j = 0; j < dimTrait; ++j) {
            displacements[pio + j] =
                    optimalRates[offset + j] * diagonal1mActualizations[pio + j];
        }
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
                double qi  = diagonal1mActualizations[offsetActualization + i];
                double qj  = diagonal1mActualizations[offsetActualization + j];

                if (Double.isInfinite(var) || (qi + qj == 0.0)) {
                    destination[destinationOffset + i * dim + j] =
                            edgeLength * variance[varianceOffset + i * dim + j];
                } else {
                    destination[destinationOffset + i * dim + j] =
                            var * (-qi * qj + qi + qj);
                }
            }
        }
    }

    /* **********************************************************************
     * Tree traversal – diagonal case
     * ******************************************************************* */

    @Override
    void actualizePrecision(DenseMatrix64F Pjp, DenseMatrix64F QjPjp,
                            int jbo, int jmo, int jdo) {

        System.arraycopy(diagonal1mActualizations,
                jdo,
                vectorDiagQdj, 0,
                dimTrait);
        oneMinus(vectorDiagQdj);

        MissingOps.diagMult(vectorDiagQdj, Pjp, QjPjp);
        MissingOps.diagMult(QjPjp,       vectorDiagQdj, Pjp);
    }

    @Override
    void actualizeVariance(DenseMatrix64F Vip, int ibo, int imo, int ido) {

        System.arraycopy(diagonal1mActualizations,
                ido,
                vectorDiagQdi, 0,
                dimTrait);
        oneMinus(vectorDiagQdi);

        diagonalDoubleProduct(Vip, vectorDiagQdi, Vip);
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] =
                    (1.0 - diagonal1mActualizations[ido + g]) * preOrderPartials[ibo + g]
                            + displacements[ido + g];
        }
    }

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip,
                                 DenseMatrix64F Pjp,
                                 DenseMatrix64F Pk) {

        System.arraycopy(diagonal1mActualizations, ido, vectorDiagQdi, 0, dimTrait);
        System.arraycopy(diagonal1mActualizations, jdo, vectorDiagQdj, 0, dimTrait);

        oneMinus(vectorDiagQdi);
        oneMinus(vectorDiagQdj);

        final DenseMatrix64F QdiPipQdi = matrix0;
        final DenseMatrix64F QdjPjpQdj = matrix1;

        diagonalDoubleProduct(Pip, vectorDiagQdi, QdiPipQdi);
        diagonalDoubleProduct(Pjp, vectorDiagQdj, QdjPjpQdj);
        CommonOps.add(QdiPipQdi, QdjPjpQdj, Pk);

        if (DEBUG) {
            System.err.println("Qdi: " + Arrays.toString(vectorDiagQdi));
            System.err.println("\tQdiPipQdi: " + QdiPipQdi);
            System.err.println("\tQdj: " + Arrays.toString(vectorDiagQdj));
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

    private static void diagonalDoubleProduct(DenseMatrix64F source,
                                              double[] diagonalScales,
                                              DenseMatrix64F destination) {
        MissingOps.diagMult(source,         diagonalScales, destination);
        MissingOps.diagMult(diagonalScales, destination);
    }
}
