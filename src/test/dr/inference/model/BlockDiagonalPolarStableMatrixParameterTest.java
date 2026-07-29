package test.dr.inference.model;

import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class BlockDiagonalPolarStableMatrixParameterTest extends MathTestCase {

    public void testAscendingOrderedRhoUsesPositiveCumulativeIncrements() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        BlockDiagonalPolarStableMatrixParameter.RhoOrdering.ASCENDING);

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
        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        BlockDiagonalPolarStableMatrixParameter.RhoOrdering.DESCENDING);

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

    public void testDefaultRhoOrderingUsesDirectRhoValues() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        BlockDiagonalPolarStableMatrixParameter.RhoOrdering.NONE);

        final double[] blockEntries = new double[blockParameter.getTridiagonalDDimension()];
        blockParameter.fillBlockDiagonalElements(blockEntries);

        assertEquals(0.50, blockEntries[0], 1e-15);
        assertEquals(0.50, blockEntries[1], 1e-15);
        assertEquals(0.25, blockEntries[2], 1e-15);
        assertEquals(0.25, blockEntries[3], 1e-15);

        final double[] nativeGradient = chainDiagonalGradient(blockParameter);
        assertEquals(3.0, nativeGradient[0], 1e-15);
        assertEquals(7.0, nativeGradient[1], 1e-15);
    }

    private static BlockDiagonalPolarStableMatrixParameter makeFourDimensionalBlockParameter(
            final double[] rho,
            final BlockDiagonalPolarStableMatrixParameter.RhoOrdering rhoOrdering) {

        final MatrixParameter rotation = makeMatrixParameter("rotation", new double[][]{
                {1.0, 0.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0},
                {0.0, 0.0, 0.0, 1.0}
        });
        final Parameter scalar = new Parameter.Default(0);
        final Parameter theta = new Parameter.Default("theta", new double[]{0.0, 0.0});
        final Parameter t = new Parameter.Default("t", new double[]{0.0, 0.0});

        return new BlockDiagonalPolarStableMatrixParameter(
                "block", rotation, scalar, new Parameter.Default("rho", rho), theta, t, rhoOrdering);
    }

    private static MatrixParameter makeMatrixParameter(final String name,
                                                       final double[][] values) {
        final int rows = values.length;
        final int cols = values[0].length;
        final Parameter[] columns = new Parameter[cols];
        for (int j = 0; j < cols; j++) {
            final double[] column = new double[rows];
            for (int i = 0; i < rows; i++) {
                column[i] = values[i][j];
            }
            columns[j] = new Parameter.Default(name + ".col" + j, column);
        }
        return new MatrixParameter(name, columns);
    }

    private static double[] chainDiagonalGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter) {

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
