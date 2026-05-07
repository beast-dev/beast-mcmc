package dr.inference.model;

import dr.xml.*;

import java.util.Collections;

/**
 * Orthogonal matrix parameterized by a sequence of Givens rotation angles.
 *
 * <p>The matrix is Q = G_0 · G_1 · ... · G_{n-1}, where each G_k is a Givens rotation
 * in the (pairI[k], pairJ[k]) coordinate plane with angle θ_k:
 * <pre>
 *   G_k[i,i] =  cos θ_k,  G_k[i,j] = -sin θ_k
 *   G_k[j,i] =  sin θ_k,  G_k[j,j] =  cos θ_k
 * </pre>
 * </p>
 */
public final class GivensRotationMatrixParameter extends AbstractComputedCompoundMatrix
        implements MatrixParameterInterface, OrthogonalMatrixProvider {

    public static final String NAME = "givensRotationMatrixParameter";

    private static final String DIMENSION = "dimension";
    private static final String ANGLES = "angles";

    private final Parameter angleParameter;
    private final int[] pairI;
    private final int[] pairJ;
    private final double[] matrixData;
    private final double[] transposeData;
    // Scratch for the 2-D → flat conversion in pullBackGradient
    private final double[] gradientMatrixScratch;
    // Rolling prefix and suffix for gradient computation — each dim×dim (O(d²) total)
    private final double[] rollingPrefixScratch;
    private final double[] rollingSuffixScratch;
    // Temp vectors used inside the per-angle gradient step — each length dim
    private final double[] tempIScratch;
    private final double[] tempJScratch;

    public GivensRotationMatrixParameter(final String name,
                                         final int dimension,
                                         final Parameter angleParameter) {
        super(name, dimension, Collections.singletonList(angleParameter));
        final int expectedAngles = dimension * (dimension - 1) / 2;
        if (angleParameter.getDimension() != expectedAngles) {
            throw new IllegalArgumentException(
                    "Expected " + expectedAngles + " Givens angles for dimension " + dimension
                            + " but found " + angleParameter.getDimension());
        }
        this.angleParameter = angleParameter;
        this.pairI = new int[expectedAngles];
        this.pairJ = new int[expectedAngles];
        this.matrixData = new double[dimension * dimension];
        this.transposeData = new double[dimension * dimension];
        this.gradientMatrixScratch = new double[dimension * dimension];
        this.rollingPrefixScratch = new double[dimension * dimension];
        this.rollingSuffixScratch = new double[dimension * dimension];
        this.tempIScratch = new double[dimension];
        this.tempJScratch = new double[dimension];

        int index = 0;
        for (int i = 0; i < dimension - 1; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                pairI[index] = i;
                pairJ[index] = j;
                index++;
            }
        }
    }

    // -------------------------------------------------------------------------
    // AbstractComputedCompoundMatrix
    // -------------------------------------------------------------------------

    @Override
    protected double computeEntry(final int row, final int col) {
        return matrixData[row * dim + col];
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return matrixData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        fillIdentity(matrixData, dim);

        for (int k = 0; k < angleParameter.getDimension(); ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);

            for (int row = 0; row < dim; ++row) {
                final int base = row * dim;
                final double oldI = matrixData[base + i];
                final double oldJ = matrixData[base + j];
                matrixData[base + i] =  oldI * c + oldJ * s;
                matrixData[base + j] = -oldI * s + oldJ * c;
            }
        }

        transpose(matrixData, transposeData, dim);
    }

    // -------------------------------------------------------------------------
    // OrthogonalMatrixProvider
    // -------------------------------------------------------------------------

    @Override
    public Parameter getOrthogonalParameter() {
        return angleParameter;
    }

    @Override
    public void fillOrthogonalMatrix(final double[] outRowMajor) {
        getParameterValue(0, 0);
        System.arraycopy(matrixData, 0, outRowMajor, 0, matrixData.length);
    }

    @Override
    public void fillOrthogonalTranspose(final double[] outRowMajor) {
        getParameterValue(0, 0);
        System.arraycopy(transposeData, 0, outRowMajor, 0, transposeData.length);
    }

    @Override
    public void applyOrthogonal(final double[] in, final double[] out) {
        System.arraycopy(in, 0, out, 0, dim);
        for (int k = angleParameter.getDimension() - 1; k >= 0; --k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);
            final double xi = out[i];
            final double xj = out[j];
            out[i] = c * xi - s * xj;
            out[j] = s * xi + c * xj;
        }
    }

    @Override
    public void applyOrthogonalTranspose(final double[] in, final double[] out) {
        System.arraycopy(in, 0, out, 0, dim);
        for (int k = 0; k < angleParameter.getDimension(); ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);
            final double xi = out[i];
            final double xj = out[j];
            out[i] =  c * xi + s * xj;
            out[j] = -s * xi + c * xj;
        }
    }

    // -------------------------------------------------------------------------
    // Gradient pull-back
    // -------------------------------------------------------------------------

    @Override
    public synchronized double[] pullBackGradient(final double[][] gradientWrtMatrix) {
        final int angleCount = angleParameter.getDimension();
        final double[] gradients = new double[angleCount];
        for (int i = 0; i < dim; ++i) {
            if (gradientWrtMatrix[i] == null || gradientWrtMatrix[i].length != dim) {
                throw new IllegalArgumentException(
                        "gradientWrtMatrix must be square with dimension " + dim);
            }
            System.arraycopy(gradientWrtMatrix[i], 0, gradientMatrixScratch, i * dim, dim);
        }
        fillPullBackGradientFlat(gradientMatrixScratch, dim, gradients);
        return gradients;
    }

    @Override
    public synchronized void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                                      final int dimension,
                                                      final double[] out,
                                                      final int offset) {
        if (dimension != dim
                || gradientWrtMatrixRowMajor == null
                || gradientWrtMatrixRowMajor.length != dim * dim) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dim * dim));
        }
        final int angleCount = angleParameter.getDimension();
        if (out == null || out.length < offset + angleCount) {
            throw new IllegalArgumentException(
                    "out must have room for " + angleCount + " entries at offset " + offset);
        }
        fillPullBackGradientRolling(gradientWrtMatrixRowMajor, out, offset);
    }

    /**
     * Computes dL/dθ_k for every Givens angle k by maintaining rolling prefix and
     * suffix matrices, each updated in O(d) per step.
     *
     * <p>Mathematical derivation:
     * <pre>
     *   Q = G_0 · ... · G_{n-1}
     *   dQ/dθ_k = prefix_k · (dG_k/dθ_k) · suffix_{k+1}
     *
     *   dL/dθ_k = Σ_{a,b} gradient[a,b] · (dQ/dθ_k)[a,b]
     *           = Σ_{c,d} (dG_k/dθ_k)[c,d] · C_k[c,d]
     *
     *   where  C_k = prefix_k^T · gradient · suffix_{k+1}^T
     * </pre>
     *
     * Because dG_k/dθ_k is zero everywhere except
     * (i,i)=−s, (i,j)=−c, (j,i)=+c, (j,j)=−s, only four entries of C_k matter.
     * They are computed by projecting onto columns i and j of the rolling prefix,
     * then dotting with rows i and j of the rolling suffix — all in O(d) after
     * a single O(d²) pass over the gradient.</p>
     *
     * <p>Memory: two dim×dim rolling matrices (O(d²)) instead of the full prefix/suffix
     * chain (O(d⁴) for O(d²) angles).</p>
     */
    private void fillPullBackGradientRolling(final double[] gradient,
                                             final double[] out,
                                             final int offset) {
        getParameterValue(0, 0); // ensure matrixData is current

        final int angleCount = angleParameter.getDimension();

        // prefix_0 = I,  suffix_0 = Q = G_0 · ... · G_{n-1}
        fillIdentity(rollingPrefixScratch, dim);
        System.arraycopy(matrixData, 0, rollingSuffixScratch, 0, dim * dim);

        for (int k = 0; k < angleCount; ++k) {
            final int pi    = pairI[k];
            final int pj    = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c  = Math.cos(theta);
            final double s  = Math.sin(theta);

            // ---- Step 1: advance suffix ----------------------------------------
            // suffix_{k+1} = G_k^T · suffix_k
            // G_k^T acts on rows pi and pj only → O(d) update.
            for (int col = 0; col < dim; ++col) {
                final double oldI = rollingSuffixScratch[pi * dim + col];
                final double oldJ = rollingSuffixScratch[pj * dim + col];
                rollingSuffixScratch[pi * dim + col] =  c * oldI + s * oldJ;
                rollingSuffixScratch[pj * dim + col] = -s * oldI + c * oldJ;
            }

            // ---- Step 2: project gradient onto columns pi and pj of prefix -----
            // tempI[q] = (prefix_k^T · gradient)[pi, q]
            //          = Σ_p  prefix_k[p, pi] · gradient[p, q]
            // tempJ[q] = same for column pj
            // Cost: O(d²) — unavoidable since gradient is dense.
            for (int q = 0; q < dim; ++q) {
                double sumI = 0.0;
                double sumJ = 0.0;
                for (int p = 0; p < dim; ++p) {
                    final double g = gradient[p * dim + q];
                    sumI += rollingPrefixScratch[p * dim + pi] * g;
                    sumJ += rollingPrefixScratch[p * dim + pj] * g;
                }
                tempIScratch[q] = sumI;
                tempJScratch[q] = sumJ;
            }

            // ---- Step 3: compute the four needed entries of C_k ----------------
            // C_k[row, col] = Σ_q  tempRow[q] · suffix_{k+1}[col, q]
            // Only rows pi and pj of suffix are accessed → O(d) shared loop.
            double Cii = 0.0, Cij = 0.0, Cji = 0.0, Cjj = 0.0;
            final int rowI = pi * dim;
            final int rowJ = pj * dim;
            for (int q = 0; q < dim; ++q) {
                final double ti = tempIScratch[q];
                final double tj = tempJScratch[q];
                final double si = rollingSuffixScratch[rowI + q];
                final double sj = rollingSuffixScratch[rowJ + q];
                Cii += ti * si;
                Cij += ti * sj;
                Cji += tj * si;
                Cjj += tj * sj;
            }

            // ---- Step 4: contract with dG_k/dθ_k ------------------------------
            // Non-zero entries: (i,i)=−s, (i,j)=−c, (j,i)=+c, (j,j)=−s
            out[offset + k] = -s * Cii - c * Cij + c * Cji - s * Cjj;

            // ---- Step 5: advance prefix ----------------------------------------
            // prefix_{k+1} = prefix_k · G_k
            // G_k acts on columns pi and pj only → O(d) update.
            for (int row = 0; row < dim; ++row) {
                final int base = row * dim;
                final double oldI = rollingPrefixScratch[base + pi];
                final double oldJ = rollingPrefixScratch[base + pj];
                rollingPrefixScratch[base + pi] =  c * oldI + s * oldJ;
                rollingPrefixScratch[base + pj] = -s * oldI + c * oldJ;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Misc public API
    // -------------------------------------------------------------------------

    public Parameter getAngleParameter() {
        return angleParameter;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    private static void fillIdentity(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension * dimension; ++i) {
            matrix[i] = 0.0;
        }
        for (int i = 0; i < dimension; ++i) {
            matrix[i * dimension + i] = 1.0;
        }
    }

    private static void transpose(final double[] in, final double[] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i * dimension + j] = in[j * dimension + i];
            }
        }
    }

    // -------------------------------------------------------------------------
    // XML parsing
    // -------------------------------------------------------------------------

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            final int dimension = xo.getIntegerAttribute(DIMENSION);
            final Parameter angles = (Parameter) xo.getElementFirstChild(ANGLES);
            return new GivensRotationMatrixParameter(name, dimension, angles);
        }

        @Override
        public String getParserDescription() {
            return "An orthogonal matrix parameterized by Givens angles.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return GivensRotationMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(ANGLES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}