package dr.inference.model;

import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Block-diagonal cos/sin parametrisation of a matrix
 *
 *      A = R D R^{-1},
 *
 * where D consists of optional 1×1 and 2×2 real blocks:
 *
 *      [ a                 r (cosθ − sinθ) ]
 *      [ r (cosθ + sinθ)             a     ]
 *
 * INTERNAL REPRESENTATION POLICY:
 *  - R, R^{-1}, A are stored as double[] (row-major)
 *  - DenseMatrix64F.wrap is used for inversion (no copying)
 *  - D is NEVER materialized densely
 *
 * @author Filippo Monti
 */
public final class BlockDiagonalCosSinMatrixParameter
        extends AbstractComputedCompoundMatrix
        implements MatrixParameterInterface {

    // ------------------------------------------------------------------
    // Dimensions and structure
    // ------------------------------------------------------------------

    private final int numBlocks;
    private final int num2x2Blocks;
    private final int[] blockStarts;
    private final int[] blockSizes;

    // ------------------------------------------------------------------
    // Underlying parameters
    // ------------------------------------------------------------------

    private final MatrixParameter RParam;
    private final Parameter aParam;
    private final Parameter rParam;
    private final Parameter thetaParam;

    // ------------------------------------------------------------------
    // Cached numeric data (authoritative)
    // ------------------------------------------------------------------

    private final double[] rData;        // R (row-major)
    private final double[] rinvData;     // R^{-1}
    private final double[] aData;        // A = R D R^{-1}

    private final double[] savedAData;
    private final double[] savedRData;
    private final double[] savedRinvData;

    private boolean rAndRinvKnown = false;
    private boolean compositionKnown = false;


    // ------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------

    public BlockDiagonalCosSinMatrixParameter(String name,
                                              MatrixParameter RParam,
                                              Parameter aParam,
                                              Parameter rParam,
                                              Parameter thetaParam) {
        super(name,
                RParam.getRowDimension(),
                Arrays.asList(RParam, aParam, rParam, thetaParam));

        this.RParam = RParam;
        this.aParam = aParam;
        this.rParam = rParam;
        this.thetaParam = thetaParam;

        if (RParam.getRowDimension() != RParam.getColumnDimension()) {
            throw new IllegalArgumentException("R must be square");
        }

        if (dim % 2 == 0) {
            num2x2Blocks = dim / 2;
            numBlocks = num2x2Blocks;
        } else {
            num2x2Blocks = (dim - 1) / 2;
            numBlocks = 1 + num2x2Blocks;
        }

        if (aParam.getDimension() != numBlocks)
            throw new IllegalArgumentException("Invalid aParam length");
        if (rParam.getDimension() != num2x2Blocks)
            throw new IllegalArgumentException("Invalid rParam length");
        if (thetaParam.getDimension() != num2x2Blocks)
            throw new IllegalArgumentException("Invalid thetaParam length");

        // Block metadata
        blockStarts = new int[numBlocks];
        blockSizes  = new int[numBlocks];
        initBlockStructure();

        // Numeric storage
        rData    = new double[dim * dim];
        rinvData = new double[dim * dim];
        aData    = new double[dim * dim];

        savedAData    = new double[dim * dim];
        savedRData    = new double[dim * dim];
        savedRinvData = new double[dim * dim];

    }

    // ------------------------------------------------------------------
    // Block structure
    // ------------------------------------------------------------------

    private void initBlockStructure() {
        if (dim % 2 == 0) {
            for (int b = 0; b < numBlocks; b++) {
                blockStarts[b] = 2 * b;
                blockSizes[b]  = 2;
            }
        } else {
            blockStarts[0] = 0;
            blockSizes[0]  = 1;
            for (int b = 0; b < num2x2Blocks; b++) {
                blockStarts[b + 1] = 1 + 2 * b;
                blockSizes[b + 1]  = 2;
            }
        }
    }

    // ------------------------------------------------------------------
    // AbstractComputedCompoundMatrix implementation
    // ------------------------------------------------------------------

    @Override
    protected double computeEntry(int row, int col) {
        ensureRAndRinvUpToDate();
        ensureCompositionUpToDate();
        return aData[row * dim + col];
    }

    @Override
    protected double getCachedValue(int row, int col) {
        return aData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        ensureRAndRinvUpToDate();
        ensureCompositionUpToDate();
    }

    // ------------------------------------------------------------------
    // Cache maintenance
    // ------------------------------------------------------------------

    private void ensureRAndRinvUpToDate() {
        if (rAndRinvKnown) return;
//        if (rAndRinvKnown) {
//            // DEBUG: detect silent changes in RParam
//            double probe = RParam.getParameterValue(0,0);
//            if (probe != rData[0]) {
//                System.err.println("WARNING: RParam changed but rAndRinvKnown is still true (missing event).");
//            }
//            return;
//        }

        // Recompute from RParam
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                rData[i * dim + j] = RParam.getParameterValue(i, j);
            }
        }

        DenseMatrix64F Rmat = DenseMatrix64F.wrap(dim, dim, rData);
        DenseMatrix64F RinvMat = DenseMatrix64F.wrap(dim, dim, rinvData);
        CommonOps.invert(Rmat, RinvMat);

        rAndRinvKnown = false; // TODO this should be true; however, rAndRinvKnown is not turned off correctly when Rparam changes (some silent stuff happening)
    }

    private void ensureCompositionUpToDate() {
        if (compositionKnown) return;

        computeA();
        compositionKnown = true;
    }

    private void computeA() {
        ensureRAndRinvUpToDate();

        final double[] aVals     = aParam.getParameterValues();
        final double[] rVals     = rParam.getParameterValues();
        final double[] thetaVals = thetaParam.getParameterValues();

        // temp = D * R^{-1}
        double[] temp = new double[dim * dim];

        int idxA = 0, idx2 = 0;

        for (int b = 0; b < numBlocks; b++) {
            int s = blockStarts[b];
            int size = blockSizes[b];

            if (size == 1) {
                double a = aVals[idxA++];
                for (int j = 0; j < dim; j++) {
                    temp[s * dim + j] += a * rinvData[s * dim + j];
                }
            } else {
                double a = aVals[idxA++];
                double r = rVals[idx2];
                double t = thetaVals[idx2++];
                double c = Math.cos(t), snt = Math.sin(t);

                double u = r * (c - snt);
                double l = r * (c + snt);

                int i0 = s, i1 = s + 1;

                for (int j = 0; j < dim; j++) {
                    double r0 = rinvData[i0 * dim + j];
                    double r1 = rinvData[i1 * dim + j];

                    temp[i0 * dim + j] += a * r0 + u * r1;
                    temp[i1 * dim + j] += l * r0 + a * r1;
                }
            }
        }

        // A = R * temp
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double sum = 0.0;
                for (int k = 0; k < dim; k++) {
                    sum += rData[i * dim + k] * temp[k * dim + j];
                }
                aData[i * dim + j] = sum;
            }
        }
    }

    // ------------------------------------------------------------------
    // Event handling override
    // ------------------------------------------------------------------

    @Override
    public void fireParameterChangedEvent() {
        rAndRinvKnown = false;
        compositionKnown = false;
        super.fireParameterChangedEvent();
    }

    // ------------------------------------------------------------------
    // Public API
    // ------------------------------------------------------------------

    public void fillRAndRinv(double[] outR, double[] outRinv) {
        ensureRAndRinvUpToDate();
        System.arraycopy(rData, 0, outR, 0, dim * dim);
        System.arraycopy(rinvData, 0, outRinv, 0, dim * dim);
    }

    public void fillBlockDiagonalElements(double[] out) {
        int expected = getCompressedDDimension();
        if (out.length != expected)
            throw new IllegalArgumentException("Invalid output length");

        for (int i = 0; i < expected; i++) out[i] = 0.0;

        final int upper = dim;
        final int lower = dim + dim - 1;

        double[] aVals = aParam.getParameterValues();
        double[] rVals = rParam.getParameterValues();
        double[] tVals = thetaParam.getParameterValues();

        int ia = 0, ib = 0;

        for (int b = 0; b < numBlocks; b++) {
            int s = blockStarts[b];
            int size = blockSizes[b];

            if (size == 1) {
                out[s] = aVals[ia++];
            } else {
                double a = aVals[ia++];
                double r = rVals[ib];
                double t = tVals[ib++];
                double c = Math.cos(t), snt = Math.sin(t);

                out[s]     = a;
                out[s + 1] = a;
                out[upper + s] = r * (c - snt);
                out[lower + s] = r * (c + snt);
            }
        }
    }

    public int[] getBlockStarts() {return blockStarts;}
    public int[] getBlockSizes() {return blockSizes;}
    public int getNumBlocks() {return numBlocks;}
    public int getNum2x2Blocks() {return num2x2Blocks;}

    public int getCompressedDDimension() {return dim + 2 * num2x2Blocks;}

    // ------------------------------------------------------------------
    // MatrixParameterInterface - delegated to parent where possible
    // ------------------------------------------------------------------

//    @Override
//    public double getParameterValue(int row, int col) {
//        // Delegated to parent class
//        return super.getParameterValue(row, col);
//    }

    public CompoundParameter getParameter() {
        CompoundParameter compoundParameter = new CompoundParameter("BlockDiagonalCosSinMatrixParameter");
        compoundParameter.addParameter(aParam);
        compoundParameter.addParameter(rParam);
        compoundParameter.addParameter(thetaParam);
        compoundParameter.addParameter(RParam);
        return compoundParameter;
    }

    public int getNumberOfParameters() {
        return aParam.getDimension() + rParam.getDimension() + thetaParam.getDimension() + RParam.getUniqueParameterCount();
    }

//    @Override
//    public double[][] getParameterAsMatrix() {
//        // Delegated to parent class
//        return super.getParameterAsMatrix();
//    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }

    // ------------------------------------------------------------------
    // MCMC state management
    // ------------------------------------------------------------------

    @Override
    protected void storeValues() {
        System.err.println("!!! BlockDiagonalCosSinMatrixParameter.storeValues() CALLED !!!");
        super.storeValues();

        if (!compositionKnown) {
            ensureRAndRinvUpToDate();
            ensureCompositionUpToDate();
        }

//        System.arraycopy(aData, 0, savedAData, 0, dim * dim);
//        System.arraycopy(rData, 0, savedRData, 0, dim * dim);
        System.arraycopy(rinvData, 0, savedRinvData, 0, dim * dim);

//        savedRAndRinvKnown = rAndRinvKnown;
//        savedCompositionKnown = compositionKnown;
    }

    @Override
    protected void restoreValues() {
        System.err.println("!!! BlockDiagonalCosSinMatrixParameter.restoreValues() CALLED !!!");
        System.err.println("    rAndRinvKnown before = " + rAndRinvKnown);

        super.restoreValues();

//        System.arraycopy(savedAData, 0, aData, 0, dim * dim);
//        System.arraycopy(savedRData, 0, rData, 0, dim * dim);
        System.arraycopy(savedRinvData, 0, rinvData, 0, dim * dim);

        rAndRinvKnown = false;
        compositionKnown = false;

        System.err.println("    rAndRinvKnown after = " + rAndRinvKnown);
    }

    /**
     * Chain rule from gradients wrt D (diagonal + first super/sub diagonals)
     * to gradients wrt (a, r, theta).
     *
     * source layout:
     *   [ 0 .. dim-1 ]          : dL/d D_{ii}
     *   [ dim .. dim+(dim-1)-1 ]: dL/d D_{i,i+1}  (upper diagonal)
     *   [ ... ]                 : dL/d D_{i+1,i}  (lower diagonal)
     *
     * out layout:
     *   [ 0 .. numBlocks-1 ]                    : dL/d a_block
     *   [ numBlocks .. numBlocks+num2x2Blocks-1 ]          : dL/d r_block
     *   [ numBlocks+num2x2Blocks .. end ]                  : dL/d theta_block
     */
    public void chainGradient(double[] source, double[] out) {

        int expectedSource = dim + 2 * num2x2Blocks;
        if (source.length != expectedSource) {
            throw new IllegalArgumentException(
                    "BlockDiagonalCosSinMatrixParameter.chainGradient: unexpected source length " +
                            source.length + " (expected " + expectedSource + ")");
        }

        int expectedOut = numBlocks + 2 * num2x2Blocks;
        if (out.length != expectedOut) {
            throw new IllegalArgumentException(
                    "BlockDiagonalCosSinMatrixParameter.chainGradient: unexpected out length " +
                            out.length + " (expected " + expectedOut + ")");
        }

        // Extract views from source (sparse format)
        double[] gDiag = source;  // first dim entries

        // Current parameter values for r and theta
        double[] rVals     = rParam.getParameterValues();
        double[] thetaVals = thetaParam.getParameterValues();

        int idx2x2 = 0; // index within 2×2 blocks

        // Layout in out[]
        int baseA     = 0;
        int baseR     = baseA + numBlocks;
        int baseTheta = baseR + num2x2Blocks;
        final int baseUpper = dim;
        final int baseLower = dim + num2x2Blocks;

        for (int b = 0; b < numBlocks; ++b) {
            int start = blockStarts[b];
            int size  = blockSizes[b];

            // --- a-block gradient (works for both 1×1 and 2×2) ---
            out[baseA + b] = chainGradientABlock(start, size, gDiag);

            // --- r, theta only for 2×2 blocks ---
            if (size == 2) {
                double r     = rVals[idx2x2];
                double theta = thetaVals[idx2x2];

                double gU = source[baseUpper + idx2x2];
                double gL = source[baseLower + idx2x2];

                out[baseR + idx2x2]     = chainGradientRBlock(gU, gL, theta);
                out[baseTheta + idx2x2] = chainGradientThetaBlock(gU, gL, r, theta);

                idx2x2++;
            }
        }
    }

    /**
     * dL/da for a single block.
     *  - 1×1 block: D_{s,s} = a      → dL/da = gDiag[s]
     *  - 2×2 block: D_{s,s} = D_{s+1,s+1} = a → dL/da = gDiag[s] + gDiag[s+1]
     */
    private double chainGradientABlock(int start, int size, double[] gDiag) {
        if (size == 1) {
            return gDiag[start];
        } else if (size == 2) {
            return gDiag[start] + gDiag[start + 1];
        } else {
            throw new IllegalArgumentException("Unsupported block size: " + size);
        }
    }

    /**
     * dL/dr for a 2×2 block.
     * Now takes the gradient values directly (already extracted from sparse array).
     */
    private double chainGradientRBlock(double gUpper, double gLower, double theta) {
        double cos = Math.cos(theta);
        double sin = Math.sin(theta);
        return gUpper * (cos - sin) + gLower * (cos + sin);
    }

    /**
     * dL/dθ for a 2×2 block.
     * Now takes the gradient values directly (already extracted from sparse array).
     */
    private double chainGradientThetaBlock(double gUpper, double gLower, double r, double theta) {
        double cos = Math.cos(theta);
        double sin = Math.sin(theta);

        double dD01_dTheta = - r * (sin + cos);
        double dD10_dTheta =  r * (cos - sin);

        return gUpper * dD01_dTheta + gLower * dD10_dTheta;
    }

    // ------------------------------------------------------------------------
    // XML Parser
    // ------------------------------------------------------------------------

    public static final String NAME          = "blockDiagonalCosSinMatrixParameter";
    private static final String R_ELEMENT    = "rotationMatrix";
    private static final String A_ELEMENT    = "blockA";
    private static final String R_BLOCK_ELT  = "blockR";
    private static final String THETA_ELEMENT= "blockTheta";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;

            MatrixParameter R =
                    (MatrixParameter) xo.getElementFirstChild(R_ELEMENT);
            Parameter aParam =
                    (Parameter) xo.getElementFirstChild(A_ELEMENT);
            Parameter rParam =
                    (Parameter) xo.getElementFirstChild(R_BLOCK_ELT);
            Parameter thetaParam =
                    (Parameter) xo.getElementFirstChild(THETA_ELEMENT);

            return new BlockDiagonalCosSinMatrixParameter(name, R, aParam, rParam, thetaParam);
        }

        @Override
        public String getParserDescription() {
            return "A matrix parameter represented as A = R D R^{-1} " +
                    "with D block-diagonal (optional 1x1 + 2x2 blocks parametrized by (a, r, theta)).";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(R_ELEMENT,     new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
                new ElementRule(A_ELEMENT,     new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(R_BLOCK_ELT,   new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(THETA_ELEMENT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        };

        @Override
        public Class getReturnType() {
            return BlockDiagonalCosSinMatrixParameter.class;
        }
    };
}