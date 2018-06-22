package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;
import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */

public class SafeMultivariateWithDriftIntegrator extends SafeMultivariateIntegrator {

    private static boolean DEBUG = false;

    public SafeMultivariateWithDriftIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
                                               int diffusionCount) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateWithDriftIntegrator");
    }

    @Override
    public void getBranchDisplacement(int bufferIndex, double[] displacement) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        System.arraycopy(displacements, bufferIndex * dimTrait,
                displacement, 0, dimTrait);
    }

    @Override
    public void getBranchExpectation(double[] actualization, double[] parentValue, double[] displacement,
                                     double[] expectation) {

        assert (expectation != null);
        assert (expectation.length >= dimTrait);

        assert (parentValue != null);
        assert (parentValue.length >= dimTrait);

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        for (int i = 0; i < dimTrait; ++i) {
            expectation[i] = parentValue[i] + displacement[i];
        }
    }

    private static final boolean TIMING = false;

    private double[] vectorDispi;
    private double[] vectorDispj;

    private void allocateStorage() {

        displacements = new double[dimTrait * bufferCount];
        vectorDispi = new double[dimTrait];
        vectorDispj = new double[dimTrait];
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Setting variances, displacement and actualization vectors
    ///////////////////////////////////////////////////////////////////////////

    public void updateBrownianDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                final double[] edgeLengths, final double[] driftRates,
                                                int updateCount) {

        super.updateBrownianDiffusionMatrices(precisionIndex, probabilityIndices, edgeLengths, driftRates, updateCount);

        if (DEBUG) {
            System.err.println("Matrices (safe with drift):");
        }

        if (driftRates != null) {

            assert (displacements != null);
            assert (driftRates.length >= updateCount * dimTrait);

            if (TIMING) {
                startTime("drift1");
            }

            int offset = 0;
            for (int up = 0; up < updateCount; ++up) {

                final double edgeLength = edgeLengths[up];
                final int pio = dimTrait * probabilityIndices[up];

                scale(driftRates, offset, edgeLength, displacements, pio, dimTrait);
                offset += dimTrait;
            }

            if (TIMING) {
                endTime("drift1");
            }
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
    void computeDelta(int jbo, int jdo, double[] delta) {
        for (int g = 0; g < dimTrait; ++g) {
            delta[g] = partials[jbo + g] - displacements[jdo + g];
        }
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] += displacements[ido + g];
        }
    }

    @Override
    InversionResult partialMean(int ibo, int jbo, int kbo,
                                int ido, int jdo) {
        if (TIMING) {
            startTime("peel4");
        }

        final double[] displacementi = vectorDispi;
        final double[] displacementj = vectorDispj;

        for (int g = 0; g < dimTrait; ++g) {
            displacementi[g] = partials[ibo + g] - displacements[ido + g];
            displacementj[g] = partials[jbo + g] - displacements[jdo + g];
        }

        final double[] tmp = vector0;

        computeWeightedSum(displacementi, displacementj, dimTrait, tmp);

        final WrappedVector kPartials = new WrappedVector.Raw(partials, kbo, dimTrait);
        final WrappedVector wrapTmp = new WrappedVector.Raw(tmp, 0, dimTrait);

        InversionResult ck = safeSolve(matrixPk, wrapTmp, kPartials, true);

        if (TIMING) {
            endTime("peel4");
            startTime("peel5");
        }
        if (DEBUG) {
            System.err.print("\t\tdisp i:");
            for (int e = 0; e < dimTrait; ++e) {
                System.err.print(" " + displacements[ido + e]);
            }
            System.err.println("");
            System.err.print("\t\tdisp j:");
            for (int e = 0; e < dimTrait; ++e) {
                System.err.print(" " + displacements[jdo + e]);
            }
        }
        return ck;
    }

    @Override
    double computeSS(final int ibo,
                     final DenseMatrix64F Pip,
                     final int jbo,
                     final DenseMatrix64F Pjp,
                     final int kbo,
                     final DenseMatrix64F Pk,
                     final int dimTrait) {
        return weightedThreeInnerProduct(vectorDispi, 0, Pip,
                vectorDispj, 0, Pjp,
                partials, kbo, Pk,
                dimTrait);
    }

    void computeWeightedSum(final double[] ipartial,
                            final double[] jpartial,
                            final int dimTrait,
                            final double[] out) {
        weightedSum(ipartial, 0, matrixPip, jpartial, 0, matrixPjp, dimTrait, out);
    }

    double[] displacements;
}
