package test.dr.inference.model;

import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class OrthogonalBlockDiagonalPolarStableMatrixParameterTest extends MathTestCase {

    public void testAscendingOrderedRhoUsesPositiveCumulativeIncrements() {
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        OrthogonalBlockDiagonalPolarStableMatrixParameter.RhoOrdering.ASCENDING);

        final double[] blockEntries = new double[blockParameter.getTridiagonalDDimension()];
        blockParameter.fillBlockDiagonalElements(blockEntries);

        assertEquals(0.50, blockEntries[0], 1e-15);
        assertEquals(0.50, blockEntries[1], 1e-15);
        assertEquals(0.75, blockEntries[2], 1e-15);
        assertEquals(0.75, blockEntries[3], 1e-15);

        final double[] nativeGradient = chainDiagonalGradient(blockParameter);
        assertEquals(10.0, nativeGradient[0], 1e-15);
        assertEquals(7.0, nativeGradient[1], 1e-15);
    }

    public void testDescendingOrderedRhoUsesPositiveCumulativeIncrements() {
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        OrthogonalBlockDiagonalPolarStableMatrixParameter.RhoOrdering.DESCENDING);

        final double[] blockEntries = new double[blockParameter.getTridiagonalDDimension()];
        blockParameter.fillBlockDiagonalElements(blockEntries);

        assertEquals(0.75, blockEntries[0], 1e-15);
        assertEquals(0.75, blockEntries[1], 1e-15);
        assertEquals(0.25, blockEntries[2], 1e-15);
        assertEquals(0.25, blockEntries[3], 1e-15);

        final double[] nativeGradient = chainDiagonalGradient(blockParameter);
        assertEquals(3.0, nativeGradient[0], 1e-15);
        assertEquals(10.0, nativeGradient[1], 1e-15);
    }

    private static OrthogonalBlockDiagonalPolarStableMatrixParameter makeFourDimensionalBlockParameter(
            final double[] rho,
            final OrthogonalBlockDiagonalPolarStableMatrixParameter.RhoOrdering rhoOrdering) {

        final Parameter angles = new Parameter.Default("angles", new double[6]);
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 4, angles);
        final Parameter scalar = new Parameter.Default(0);
        final Parameter theta = new Parameter.Default("theta", new double[]{0.0, 0.0});
        final Parameter t = new Parameter.Default("t", new double[]{0.0, 0.0});

        return new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                "block", rotation, scalar, new Parameter.Default("rho", rho), theta, t, rhoOrdering);
    }

    private static double[] chainDiagonalGradient(
            final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter) {

        final double[] compressedGradient = new double[blockParameter.getCompressedDDimension()];
        compressedGradient[0] = 1.0;
        compressedGradient[1] = 2.0;
        compressedGradient[2] = 3.0;
        compressedGradient[3] = 4.0;

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedGradient, nativeGradient);
        return nativeGradient;
    }
}
