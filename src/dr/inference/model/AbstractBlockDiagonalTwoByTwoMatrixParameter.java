package dr.inference.model;

import dr.math.matrixAlgebra.WrappedMatrix;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Generic block-diagonal matrix parameter of the form
 *
 *     A = R D R^{-1},
 *
 * where D has either:
 * - only 2x2 blocks   (even dimension), or
 * - one leading 1x1 block followed by 2x2 blocks (odd dimension).
 *
 * Unlike the older superclass, this class does NOT assume that a 2x2 block has the form
 *
 *     [ a  u ]
 *     [ l  a ]
 *
 * Instead, each subclass provides a mapping from native block parameters to the full raw 2x2 block
 *
 *     [ d00  d01 ]
 *     [ d10  d11 ]
 *
 * and a pullback from raw-entry gradients (g00, g01, g10, g11) to the native block parameters.
 *
 * This makes the class flexible enough to support:
 * - the old equal-diagonal charts,
 * - the Cartesian stable chart,
 * - the true polar stable chart,
 * - future 2x2 charts with unequal diagonals if needed.
 */
public abstract class AbstractBlockDiagonalTwoByTwoMatrixParameter
        extends AbstractComputedCompoundMatrix
        implements MatrixParameterInterface {

    protected final int numBlocks;
    protected final int num2x2Blocks;
    protected final boolean hasLeadingOneByOneBlock;
    protected final int[] blockStarts;
    protected final int[] blockSizes;

    protected final MatrixParameter RParam;

    /**
     * Parameter for the optional leading 1x1 block.
     * Must have dimension 1 if dim is odd, and 0 if dim is even.
     */
    protected final Parameter scalarBlockParam;

    /**
     * Native parameters for the 2x2 blocks.
     * Every parameter in this array must have dimension num2x2Blocks.
     * For example:
     * - old chart: [a, p1, p2]
     * - polar chart: [rho, theta, t]
     * - Cartesian chart: [alpha, sigma, t]
     */
    protected final Parameter[] twoByTwoBlockParams;

    protected final double[] rData;
    protected final double[] rinvData;
    protected final double[] aData;
    private final double[] temp;
    private final double[] savedRData;
    private final double[] savedRinvData;
    private boolean savedRAndRinvKnown = false;

    private boolean rAndRinvKnown = false;
    private boolean compositionKnown = false;

    private final CompoundParameter compoundParameter;

    protected AbstractBlockDiagonalTwoByTwoMatrixParameter(final String name,
                                                                  final MatrixParameter RParam,
                                                                  final Parameter scalarBlockParam,
                                                                  final Parameter... twoByTwoBlockParams) {
        super(name,
                RParam.getRowDimension(),
                buildParameterList(RParam, scalarBlockParam, twoByTwoBlockParams));

        this.RParam = RParam;
        if (scalarBlockParam != null) {
            this.scalarBlockParam = scalarBlockParam;
        } else {
            this.scalarBlockParam = null;
        }

        this.twoByTwoBlockParams = twoByTwoBlockParams.clone();

        if (RParam.getRowDimension() != RParam.getColumnDimension()) {
            throw new IllegalArgumentException("R must be square");
        }

        this.hasLeadingOneByOneBlock = (dim & 1) == 1;
        this.num2x2Blocks = dim / 2;
        this.numBlocks = hasLeadingOneByOneBlock ? num2x2Blocks + 1 : num2x2Blocks;

//        if (scalarBlockParam == null) {
//            throw new IllegalArgumentException("scalarBlockParam must not be null");
//        }

        final int expectedScalarDim = hasLeadingOneByOneBlock ? 1 : 0;
//        if (scalarBlockParam.getDimension() != expectedScalarDim) {
//            throw new IllegalArgumentException(
//                    "scalarBlockParam must have dimension " + expectedScalarDim +
//                            " but has dimension " + scalarBlockParam.getDimension());
//        }

        if (twoByTwoBlockParams.length != getTwoByTwoBlockParameterCount()) {
            throw new IllegalArgumentException(
                    "Expected " + getTwoByTwoBlockParameterCount() +
                            " two-by-two block parameters, but received " + twoByTwoBlockParams.length);
        }

        for (int k = 0; k < twoByTwoBlockParams.length; k++) {
            if (twoByTwoBlockParams[k] == null) {
                throw new IllegalArgumentException("twoByTwoBlockParams[" + k + "] is null");
            }
            if (twoByTwoBlockParams[k].getDimension() != num2x2Blocks) {
                throw new IllegalArgumentException(
                        "Invalid dimension for parameter '" + getTwoByTwoBlockParameterName(k) +
                                "': expected " + num2x2Blocks +
                                ", found " + twoByTwoBlockParams[k].getDimension());
            }
        }

        this.blockStarts = new int[numBlocks];
        this.blockSizes = new int[numBlocks];
        initBlockStructure();

        final int matrixSize = dim * dim;
        this.rData = new double[matrixSize];
        this.rinvData = new double[matrixSize];
        this.aData = new double[matrixSize];
        this.temp = new double[matrixSize];
        this.savedRData = new double[matrixSize];
        this.savedRinvData = new double[matrixSize];

        this.compoundParameter = new CompoundParameter(getClass().getSimpleName());
        this.compoundParameter.addParameter(scalarBlockParam);
        for (final Parameter p : twoByTwoBlockParams) {
            this.compoundParameter.addParameter(p);
        }
        this.compoundParameter.addParameter(RParam);
    }

    private static List<Parameter> buildParameterList(final MatrixParameter RParam,
                                                      final Parameter scalarBlockParam,
                                                      final Parameter... twoByTwoBlockParams) {
        final List<Parameter> list = new ArrayList<Parameter>(2 + twoByTwoBlockParams.length);
        list.add(RParam);
        if (scalarBlockParam != null) {
            list.add(scalarBlockParam);
        }
        list.addAll(Arrays.asList(twoByTwoBlockParams));
        return list;
    }

    private void initBlockStructure() {
        if (!hasLeadingOneByOneBlock) {
            for (int b = 0; b < numBlocks; b++) {
                blockStarts[b] = 2 * b;
                blockSizes[b] = 2;
            }
        } else {
            blockStarts[0] = 0;
            blockSizes[0] = 1;
            for (int b = 0; b < num2x2Blocks; b++) {
                blockStarts[b + 1] = 1 + 2 * b;
                blockSizes[b + 1] = 2;
            }
        }
    }

    @Override
    protected double computeEntry(final int row, final int col) {
        ensureUpToDate();
        return aData[row * dim + col];
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return aData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        ensureUpToDate();
    }

    private void ensureUpToDate() {
        ensureRAndRinvUpToDate();
        ensureCompositionUpToDate();
    }

    private void ensureRAndRinvUpToDate() {
        if (rAndRinvKnown) {
            return;
        }

        int idx = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                rData[idx++] = RParam.getParameterValue(i, j);
            }
        }

        final DenseMatrix64F Rmat = DenseMatrix64F.wrap(dim, dim, rData);
        final DenseMatrix64F RinvMat = DenseMatrix64F.wrap(dim, dim, rinvData);
        if (!CommonOps.invert(Rmat, RinvMat)) {
            throw new IllegalArgumentException("R is singular and cannot be inverted");
        }

        rAndRinvKnown = true;
    }

    private void ensureCompositionUpToDate() {
        if (compositionKnown) {
            return;
        }
        computeA();
        compositionKnown = true;
    }

    private void computeA() {
        Arrays.fill(temp, 0.0);

        final double[][] blockParamValues = getTwoByTwoBlockParameterValues();
        final double[] rawBlock = new double[4]; // [d00, d01, d10, d11]

        int block2x2Idx = 0;
        for (int b = 0; b < numBlocks; b++) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];

            if (size == 1) {
                final double s = scalarBlockParam.getParameterValue(0);
                final int rowOffset = start * dim;
                for (int j = 0; j < dim; j++) {
                    temp[rowOffset + j] = s * rinvData[rowOffset + j];
                }
                continue;
            }

            fillTwoByTwoBlock(block2x2Idx, blockParamValues, rawBlock);

            final double d00 = rawBlock[0];
            final double d01 = rawBlock[1];
            final double d10 = rawBlock[2];
            final double d11 = rawBlock[3];

            final int i0 = start;
            final int i1 = start + 1;
            final int row0 = i0 * dim;
            final int row1 = i1 * dim;

            for (int j = 0; j < dim; j++) {
                final double r0 = rinvData[row0 + j];
                final double r1 = rinvData[row1 + j];
                temp[row0 + j] = d00 * r0 + d01 * r1;
                temp[row1 + j] = d10 * r0 + d11 * r1;
            }

            block2x2Idx++;
        }

        for (int i = 0; i < dim; i++) {
            final int row = i * dim;
            for (int j = 0; j < dim; j++) {
                double sum = 0.0;
                for (int k = 0; k < dim; k++) {
                    sum += rData[row + k] * temp[k * dim + j];
                }
                aData[row + j] = sum;
            }
        }
    }

    public void fillRAndRinv(final double[] outR, final double[] outRinv) {
        ensureRAndRinvUpToDate();
        System.arraycopy(rData, 0, outR, 0, dim * dim);
        System.arraycopy(rinvData, 0, outRinv, 0, dim * dim);
    }

    /**
     * Fills the block-diagonal entries in the tridiagonal-style representation
     *
     *     [diag(0..dim-1), upper(0..dim-2), lower(0..dim-2)].
     *
     * Entries outside the fixed block structure are zero.
     */
    public void fillBlockDiagonalElements(final double[] out) {
        final int expected = getTridiagonalDDimension();
        if (out.length != expected) {
            throw new IllegalArgumentException("Invalid output length");
        }

        Arrays.fill(out, 0.0);

        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        final double[][] blockParamValues = getTwoByTwoBlockParameterValues();
        final double[] rawBlock = new double[4];

        int block2x2Idx = 0;
        for (int b = 0; b < numBlocks; b++) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];

            if (size == 1) {
                out[start] = scalarBlockParam.getParameterValue(0);
                continue;
            }

            fillTwoByTwoBlock(block2x2Idx, blockParamValues, rawBlock);

            out[start] = rawBlock[0];
            out[start + 1] = rawBlock[3];
            out[upperOffset + start] = rawBlock[1];
            out[lowerOffset + start] = rawBlock[2];

            block2x2Idx++;
        }
    }

    /**
     * Pulls back gradients from the raw block representation
     *
     * raw source layout:
     *   - diag entries: dim values
     *   - active upper entries: num2x2Blocks values
     *   - active lower entries: num2x2Blocks values
     *
     * output layout:
     *   - optional scalar 1x1 parameter first (0 or 1 value)
     *   - then all native 2x2 block parameters grouped by parameter family
     */
    public void chainGradient(final double[] source, final double[] out) {
        final int expectedSource = getCompressedDDimension();
        if (source.length != expectedSource) {
            throw new IllegalArgumentException(
                    "Unexpected source length " + source.length +
                            " (expected " + expectedSource + ")");
        }

        final int expectedOut = getBlockDiagonalNParameters();
        if (out.length != expectedOut) {
            throw new IllegalArgumentException(
                    "Unexpected out length " + out.length +
                            " (expected " + expectedOut + ")");
        }

        Arrays.fill(out, 0.0);

        final double[][] blockParamValues = getTwoByTwoBlockParameterValues();

        final int scalarBase = 0;
        final int twoByTwoBase = hasLeadingOneByOneBlock ? 1 : 0;

        final int upperBase = dim;
        final int lowerBase = dim + num2x2Blocks;

        int block2x2Idx = 0;
        for (int b = 0; b < numBlocks; b++) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];

            if (size == 1) {
                out[scalarBase] = source[start];
                continue;
            }

            final double g00 = source[start];
            final double g11 = source[start + 1];
            final double g01 = source[upperBase + block2x2Idx];
            final double g10 = source[lowerBase + block2x2Idx];

            chainTwoByTwoBlockGradient(
                    block2x2Idx,
                    g00, g01, g10, g11,
                    blockParamValues,
                    out,
                    twoByTwoBase);

            block2x2Idx++;
        }
    }

    private double[][] getTwoByTwoBlockParameterValues() {
        final int k = twoByTwoBlockParams.length;
        final double[][] values = new double[k][];
        for (int i = 0; i < k; i++) {
            values[i] = twoByTwoBlockParams[i].getParameterValues();
        }
        return values;
    }

    public int[] getBlockStarts() {
        return blockStarts;
    }

    public int[] getBlockSizes() {
        return blockSizes;
    }

    public int getNumBlocks() {
        return numBlocks;
    }

    public int getNum2x2Blocks() {
        return num2x2Blocks;
    }

    public boolean hasLeadingOneByOneBlock() {
        return hasLeadingOneByOneBlock;
    }

    /**
     * Raw source gradient dimension:
     * diag entries (dim) + active uppers + active lowers.
     */
    public int getCompressedDDimension() {
        return dim + 2 * num2x2Blocks;
    }

    /**
     * Native block parameter dimension:
     * optional leading scalar + all 2x2 native parameter families.
     */
    public int getBlockDiagonalNParameters() {
        return scalarBlockParam.getDimension() + getTwoByTwoBlockParameterCount() * num2x2Blocks;
    }

    public int getTridiagonalDDimension() {
        return 3 * dim - 2;
    }

    public int getTwoByTwoParameterFamilyCount() {
        return getTwoByTwoBlockParameterCount();
    }

    public CompoundParameter getParameter() {
        return compoundParameter;
    }

    public int getNumberOfParameters() {
        int total = scalarBlockParam.getDimension() + dim * dim;
        for (final Parameter p : twoByTwoBlockParams) {
            total += p.getDimension();
        }
        return total;
    }

    public Parameter getScalarBlockParameter() {
        return scalarBlockParam;
    }

    public Parameter getTwoByTwoBlockParameter(final int index) {
        return twoByTwoBlockParams[index];
    }

    public MatrixParameter getRotationMatrixParameter() {
        return RParam;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }

    @Override
    public void fireParameterChangedEvent() {
        rAndRinvKnown = false;
        compositionKnown = false;
        super.fireParameterChangedEvent();
    }

    @Override
    public void variableChangedEvent(final Variable variable,
                                     final int index,
                                     final Parameter.ChangeType type) {

        if (variable == RParam) {
            rAndRinvKnown = false;
        }

        if (variable == RParam || variable == scalarBlockParam) {
            compositionKnown = false;
        } else {
            for (final Parameter p : twoByTwoBlockParams) {
                if (variable == p) {
                    compositionKnown = false;
                    break;
                }
            }
        }

        super.variableChangedEvent(variable, index, type);
    }

    @Override
    protected void storeValues() {
        super.storeValues();
        if (!compositionKnown) {
            ensureRAndRinvUpToDate();
            ensureCompositionUpToDate();
        }
        System.arraycopy(rData, 0, savedRData, 0, dim * dim);
        System.arraycopy(rinvData, 0, savedRinvData, 0, dim * dim);
        savedRAndRinvKnown = rAndRinvKnown;
    }

    @Override
    protected void restoreValues() {
        super.restoreValues();
        if (savedRAndRinvKnown) {
            System.arraycopy(savedRData, 0, rData, 0, dim * dim);
            System.arraycopy(savedRinvData, 0, rinvData, 0, dim * dim);
        }
        rAndRinvKnown = savedRAndRinvKnown;
        compositionKnown = false; // A depends on restored native params, so rebuild lazily.
    }

    /**
     * Number of native parameter families used by each 2x2 block.
     * Examples:
     * - old equal-diagonal chart: 3  (a, p1, p2)
     * - polar chart: 3  (rho, theta, t)
     * - future chart with 4 params: 4
     */
    protected abstract int getTwoByTwoBlockParameterCount();

    /**
     * Human-readable name for the k-th native 2x2 block parameter family.
     */
    protected abstract String getTwoByTwoBlockParameterName(int k);

    /**
     * Build the raw 2x2 block for block index b in row-major order:
     *
     * outBlock[0] = d00
     * outBlock[1] = d01
     * outBlock[2] = d10
     * outBlock[3] = d11
     */
    protected abstract void fillTwoByTwoBlock(int blockIndex,
                                              double[][] blockParameterValues,
                                              double[] outBlock);

    /**
     * Pull back raw gradients for one 2x2 block to the native parameters.
     *
     * The out layout is grouped by parameter family:
     *
     *   out[base + 0 * num2x2Blocks + blockIndex]  <- grad wrt first block parameter
     *   out[base + 1 * num2x2Blocks + blockIndex]  <- grad wrt second block parameter
     *   ...
     *
     * where base = (hasLeadingOneByOneBlock ? 1 : 0).
     */
    protected abstract void chainTwoByTwoBlockGradient(int blockIndex,
                                                       double g00,
                                                       double g01,
                                                       double g10,
                                                       double g11,
                                                       double[][] blockParameterValues,
                                                       double[] out,
                                                       int baseOffset);

    /**
     * Helper for subclasses: index into the grouped output layout.
     */
    protected final int twoByTwoGradientOffset(final int baseOffset,
                                               final int parameterIndex,
                                               final int blockIndex) {
        return baseOffset + parameterIndex * num2x2Blocks + blockIndex;
    }
}
