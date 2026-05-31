package test.dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;

public final class BlockDiagonalLyapunovEqualDiagonalPlanTest extends TestCase {

    private static final double TOL = 1.0e-10;

    public void testAdjointHelperMatchesGenericSolverWithLeadingScalar() {
        assertAdjointMatchesGeneric(leadingScalarStructure(), leadingScalarParams(), symmetricMatrix(5));
    }

    public void testAdjointHelperMatchesGenericSolverWithoutLeadingScalar() {
        assertAdjointMatchesGeneric(allTwoByTwoStructure(), allTwoByTwoParams(), symmetricMatrix(4));
    }

    public void testAdjointHelperFallsBackForNonSymmetricRhs() {
        final BlockDiagonalExpSolver.BlockStructure structure = leadingScalarStructure();
        final double[] blockDParams = leadingScalarParams();
        final DenseMatrix64F hBlock = symmetricMatrix(structure.getDim());
        hBlock.data[1 * structure.getDim() + 4] += 0.37;

        assertAdjointMatchesGeneric(structure, blockDParams, hBlock);
    }

    public void testAdjointHelperFallsBackForUnequalBlockDiagonalParams() {
        final BlockDiagonalExpSolver.BlockStructure structure = allTwoByTwoStructure();
        final double[] blockDParams = allTwoByTwoParams();
        blockDParams[1] += 0.11;

        assertAdjointMatchesGeneric(structure, blockDParams, symmetricMatrix(structure.getDim()));
    }

    public void testAdjointHelperSupportsAliasedInputOutput() {
        final BlockDiagonalExpSolver.BlockStructure structure = allTwoByTwoStructure();
        final double[] blockDParams = allTwoByTwoParams();
        final DenseMatrix64F hBlock = symmetricMatrix(structure.getDim());

        final BlockDiagonalLyapunovAdjointHelper helper = createHelper(structure);
        final DenseMatrix64F expected = new DenseMatrix64F(structure.getDim(), structure.getDim());
        helper.solveAdjointInDBasis(hBlock, blockDParams, expected);

        final DenseMatrix64F aliased = hBlock.copy();
        helper.solveAdjointInDBasis(aliased, blockDParams, aliased);

        assertMatrixEquals(expected, aliased, TOL);
    }

    private static void assertAdjointMatchesGeneric(final BlockDiagonalExpSolver.BlockStructure structure,
                                                    final double[] blockDParams,
                                                    final DenseMatrix64F hBlock) {
        final DenseMatrix64F negH = new DenseMatrix64F(structure.getDim(), structure.getDim());
        negate(hBlock, negH);

        final DenseMatrix64F expected =
                solveGeneric(structure, transposeParams(blockDParams, structure.getDim()), negH);
        final DenseMatrix64F actual = new DenseMatrix64F(structure.getDim(), structure.getDim());

        createHelper(structure).solveAdjointInDBasis(hBlock, blockDParams, actual);

        assertMatrixEquals(expected, actual, TOL);
    }

    private static BlockDiagonalLyapunovAdjointHelper createHelper(
            final BlockDiagonalExpSolver.BlockStructure structure) {
        final BlockDiagonalLyapunovSolver solver = new BlockDiagonalLyapunovSolver(
                structure.getDim(), structure.getBlockStarts(), structure.getBlockSizes());
        return new BlockDiagonalLyapunovAdjointHelper(structure.getDim(), solver, structure);
    }

    private static BlockDiagonalExpSolver.BlockStructure leadingScalarStructure() {
        return new BlockDiagonalExpSolver.BlockStructure(
                5,
                new int[]{0, 1, 3},
                new int[]{1, 2, 2});
    }

    private static BlockDiagonalExpSolver.BlockStructure allTwoByTwoStructure() {
        return new BlockDiagonalExpSolver.BlockStructure(
                4,
                new int[]{0, 2},
                new int[]{2, 2});
    }

    private static double[] leadingScalarParams() {
        final int dim = 5;
        final double[] params = new double[3 * dim - 2];
        params[0] = 1.10;
        params[1] = 0.75;
        params[2] = 0.75;
        params[3] = 1.35;
        params[4] = 1.35;
        params[dim + 1] = 0.22;
        params[dim + 3] = -0.16;
        params[dim + (dim - 1) + 1] = -0.31;
        params[dim + (dim - 1) + 3] = 0.27;
        return params;
    }

    private static double[] allTwoByTwoParams() {
        final int dim = 4;
        final double[] params = new double[3 * dim - 2];
        params[0] = 0.80;
        params[1] = 0.80;
        params[2] = 1.20;
        params[3] = 1.20;
        params[dim] = 0.25;
        params[dim + 2] = -0.35;
        params[dim + (dim - 1)] = -0.18;
        params[dim + (dim - 1) + 2] = 0.21;
        return params;
    }

    private static DenseMatrix64F solveGeneric(final BlockDiagonalExpSolver.BlockStructure structure,
                                               final double[] blockDParams,
                                               final DenseMatrix64F rhs) {
        final DenseMatrix64F out = new DenseMatrix64F(structure.getDim(), structure.getDim());
        final BlockDiagonalLyapunovSolver solver = new BlockDiagonalLyapunovSolver(
                structure.getDim(), structure.getBlockStarts(), structure.getBlockSizes());
        solver.solveOrThrow(blockDParams, rhs, out);
        return out;
    }

    private static DenseMatrix64F symmetricMatrix(final int dim) {
        final DenseMatrix64F matrix = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = i; j < dim; ++j) {
                final double value = 0.1 * (i + 1) + 0.07 * (j + 2) + 0.03 * (i + 1) * (j + 1);
                matrix.data[i * dim + j] = value;
                matrix.data[j * dim + i] = value;
            }
        }
        return matrix;
    }

    private static double[] transposeParams(final double[] blockDParams,
                                            final int dim) {
        final int offLen = dim - 1;
        final double[] transposed = new double[blockDParams.length];
        System.arraycopy(blockDParams, 0, transposed, 0, dim);
        System.arraycopy(blockDParams, dim + offLen, transposed, dim, offLen);
        System.arraycopy(blockDParams, dim, transposed, dim + offLen, offLen);
        return transposed;
    }

    private static void negate(final DenseMatrix64F source,
                               final DenseMatrix64F out) {
        for (int i = 0; i < source.data.length; ++i) {
            out.data[i] = -source.data[i];
        }
    }

    private static void assertMatrixEquals(final DenseMatrix64F expected,
                                           final DenseMatrix64F actual,
                                           final double tolerance) {
        assertEquals(expected.numRows, actual.numRows);
        assertEquals(expected.numCols, actual.numCols);
        for (int i = 0; i < expected.data.length; ++i) {
            assertEquals("entry " + i, expected.data[i], actual.data[i], tolerance);
        }
    }
}
