package test.dr.inference.model;

import dr.inference.model.BlockDiagonalNormalizedCartesianStableMatrixParameter;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class BlockDiagonalNormalizedCartesianStableMatrixParameterTest extends MathTestCase {

    public void testBlockEntriesUsePositiveSelectionConvention() {
        final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter =
                makeTwoDimensionalBlockParameter(2.0, 0.25, 0.7);

        final double[] blockEntries = new double[blockParameter.getTridiagonalDDimension()];
        blockParameter.fillBlockDiagonalElements(blockEntries);

        assertEquals(2.0, blockEntries[0], 1e-15);
        assertEquals(2.0, blockEntries[1], 1e-15);
        assertEquals(-0.2, blockEntries[2], 1e-15);
        assertEquals(1.2, blockEntries[3], 1e-15);
    }

    public void testMatchesPolarUnderNormalizedCartesianSubstitution() {
        final double[] rate = new double[]{0.8, 1.4};
        final double[] ratio = new double[]{-0.2, 0.6};
        final double[] t = new double[]{0.3, 1.1};
        final BlockDiagonalNormalizedCartesianStableMatrixParameter cartesian =
                makeFourDimensionalBlockParameter(
                        rate,
                        ratio,
                        t,
                        BlockDiagonalNormalizedCartesianStableMatrixParameter.RateOrdering.NONE);

        final double[] rho = new double[rate.length];
        final double[] theta = new double[rate.length];
        for (int i = 0; i < rate.length; i++) {
            rho[i] = rate[i] * Math.sqrt(1.0 + ratio[i] * ratio[i]);
            theta[i] = Math.atan(ratio[i]);
        }
        final BlockDiagonalPolarStableMatrixParameter polar =
                makePolarBlockParameter(rho, theta, t);

        assertBlockEntriesEqual(polar, cartesian, 1e-15);
    }

    public void testAscendingOrderedRateUsesPositiveCumulativeIncrements() {
        final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        new double[]{0.0, 0.0},
                        new double[]{0.0, 0.0},
                        BlockDiagonalNormalizedCartesianStableMatrixParameter.RateOrdering.ASCENDING);

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

    public void testDescendingOrderedRateUsesPositiveCumulativeIncrements() {
        final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        new double[]{0.0, 0.0},
                        new double[]{0.0, 0.0},
                        BlockDiagonalNormalizedCartesianStableMatrixParameter.RateOrdering.DESCENDING);

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

    public void testDefaultRateOrderingUsesDirectRateValues() {
        final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter =
                makeFourDimensionalBlockParameter(
                        new double[]{0.50, 0.25},
                        new double[]{0.0, 0.0},
                        new double[]{0.0, 0.0},
                        BlockDiagonalNormalizedCartesianStableMatrixParameter.RateOrdering.NONE);

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

    public void testPullBackGradientMatchesFiniteDifference() {
        final Parameter rate = new Parameter.Default("rate", new double[]{0.7, 1.3});
        final Parameter ratio = new Parameter.Default("ratio", new double[]{-0.4, 0.2});
        final Parameter t = new Parameter.Default("t", new double[]{0.5, -0.8});
        final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter =
                new BlockDiagonalNormalizedCartesianStableMatrixParameter(
                        "cartesian",
                        identityMatrixParameter(4),
                        new Parameter.Default(0),
                        rate,
                        ratio,
                        t);
        final double[] compressedGradient = new double[]{
                1.1, -0.7, 0.4, 2.0,
                -1.3, 0.9,
                0.6, -0.2
        };

        final double[] analytic = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedGradient, analytic);

        final Parameter[] parameters = new Parameter[]{rate, ratio, t};
        int offset = 0;
        for (Parameter parameter : parameters) {
            for (int i = 0; i < parameter.getDimension(); i++) {
                assertEquals(
                        "finite difference gradient for " + parameter.getId() + "[" + i + "]",
                        finiteDifference(blockParameter, parameter, i, compressedGradient),
                        analytic[offset + i],
                        1e-6);
            }
            offset += parameter.getDimension();
        }
    }

    private static BlockDiagonalNormalizedCartesianStableMatrixParameter makeTwoDimensionalBlockParameter(
            final double rate,
            final double ratio,
            final double t) {

        return new BlockDiagonalNormalizedCartesianStableMatrixParameter(
                "cartesian",
                identityMatrixParameter(2),
                new Parameter.Default(0),
                new Parameter.Default("rate", new double[]{rate}),
                new Parameter.Default("ratio", new double[]{ratio}),
                new Parameter.Default("t", new double[]{t}));
    }

    private static BlockDiagonalNormalizedCartesianStableMatrixParameter makeFourDimensionalBlockParameter(
            final double[] rate,
            final double[] ratio,
            final double[] t,
            final BlockDiagonalNormalizedCartesianStableMatrixParameter.RateOrdering rateOrdering) {

        return new BlockDiagonalNormalizedCartesianStableMatrixParameter(
                "cartesian",
                identityMatrixParameter(4),
                new Parameter.Default(0),
                new Parameter.Default("rate", rate),
                new Parameter.Default("ratio", ratio),
                new Parameter.Default("t", t),
                rateOrdering);
    }

    private static BlockDiagonalPolarStableMatrixParameter makePolarBlockParameter(
            final double[] rho,
            final double[] theta,
            final double[] t) {

        return new BlockDiagonalPolarStableMatrixParameter(
                "polar",
                identityMatrixParameter(4),
                new Parameter.Default(0),
                new Parameter.Default("rho", rho),
                new Parameter.Default("theta", theta),
                new Parameter.Default("t", t));
    }

    private static MatrixParameter identityMatrixParameter(final int dim) {
        final Parameter[] columns = new Parameter[dim];
        for (int j = 0; j < dim; j++) {
            final double[] column = new double[dim];
            column[j] = 1.0;
            columns[j] = new Parameter.Default("identity.col" + j, column);
        }
        return new MatrixParameter("identity", columns);
    }

    private static void assertBlockEntriesEqual(
            final BlockDiagonalPolarStableMatrixParameter polar,
            final BlockDiagonalNormalizedCartesianStableMatrixParameter cartesian,
            final double tolerance) {

        final double[] polarEntries = new double[polar.getTridiagonalDDimension()];
        final double[] cartesianEntries = new double[cartesian.getTridiagonalDDimension()];
        polar.fillBlockDiagonalElements(polarEntries);
        cartesian.fillBlockDiagonalElements(cartesianEntries);
        assertEquals(polarEntries.length, cartesianEntries.length);
        for (int i = 0; i < polarEntries.length; i++) {
            assertEquals("block entry " + i, polarEntries[i], cartesianEntries[i], tolerance);
        }
    }

    private static double[] chainDiagonalGradient(
            final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter) {

        final double[] compressedGradient = new double[blockParameter.getCompressedDDimension()];
        compressedGradient[0] = 1.0;
        compressedGradient[1] = 2.0;
        compressedGradient[2] = 3.0;
        compressedGradient[3] = 4.0;

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedGradient, nativeGradient);
        return nativeGradient;
    }

    private static double finiteDifference(
            final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter,
            final Parameter parameter,
            final int index,
            final double[] compressedGradient) {

        final double h = 1e-7;
        final double original = parameter.getParameterValue(index);
        parameter.setParameterValue(index, original + h);
        blockParameter.fireParameterChangedEvent();
        final double plus = objective(blockParameter, compressedGradient);
        parameter.setParameterValue(index, original - h);
        blockParameter.fireParameterChangedEvent();
        final double minus = objective(blockParameter, compressedGradient);
        parameter.setParameterValue(index, original);
        blockParameter.fireParameterChangedEvent();
        return (plus - minus) / (2.0 * h);
    }

    private static double objective(
            final BlockDiagonalNormalizedCartesianStableMatrixParameter blockParameter,
            final double[] compressedGradient) {

        final double[] blockEntries = new double[blockParameter.getTridiagonalDDimension()];
        blockParameter.fillBlockDiagonalElements(blockEntries);
        double objective = 0.0;

        final int dim = blockParameter.getRowDimension();
        final int num2x2Blocks = blockParameter.getNum2x2Blocks();
        for (int i = 0; i < dim; i++) {
            objective += compressedGradient[i] * blockEntries[i];
        }

        final int compressedUpperOffset = dim;
        final int compressedLowerOffset = dim + num2x2Blocks;
        final int tridiagonalUpperOffset = dim;
        final int tridiagonalLowerOffset = dim + dim - 1;
        final int[] blockStarts = blockParameter.getBlockStarts();
        final int[] blockSizes = blockParameter.getBlockSizes();
        int block2x2Index = 0;
        for (int blockIndex = 0; blockIndex < blockStarts.length; blockIndex++) {
            if (blockSizes[blockIndex] == 1) {
                continue;
            }
            final int start = blockStarts[blockIndex];
            objective += compressedGradient[compressedUpperOffset + block2x2Index] *
                    blockEntries[tridiagonalUpperOffset + start];
            objective += compressedGradient[compressedLowerOffset + block2x2Index] *
                    blockEntries[tridiagonalLowerOffset + start];
            block2x2Index++;
        }
        return objective;
    }
}
