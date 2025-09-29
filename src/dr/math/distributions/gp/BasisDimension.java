package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;

public class BasisDimension {

    private final GaussianProcessKernel kernel;
    private final DesignMatrix design1;
    private final DesignMatrix design2;
    private final WeightFunction weightFunction;
    private double[] w1, w2;
    private DenseMatrix64F scratchK;

    public BasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2) {
        this(kernel, design1, design2, null);
    }

    public BasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2, WeightFunction weightFunction) {
        this.kernel = kernel;
        this.design1 = design1;
        this.design2 = design2;
        this.weightFunction = weightFunction;
    }

    public BasisDimension(GaussianProcessKernel kernel, DesignMatrix design) {
        this(kernel, design, design);
    }

    public BasisDimension(GaussianProcessKernel kernel, RandomField.WeightProvider weights) {
        this(kernel, makeDesignMatrixFromWeights(weights));
    }

    public final void addFirstOrderKernelComponent(DenseMatrix64F gramian) {
        final int nR = gramian.getNumRows();
        final int nC = gramian.getNumCols();
        ensureCapacity(nR, nC);

        final double[] x1 = design1.getColumnValues(0); //  TODO: generalize to multi-D row
        final double[] x2 = design2.getColumnValues(0);

        final boolean useWeights = (weightFunction != null);
        if (useWeights) {
            precomputeWeights(x1, nR, w1);
            precomputeWeights(x2, nC, w2);
        }

        // Decide target: write directly to gramian unless subclass needs a preimage K
        final boolean needsPost = needsPostProcess();
        final DenseMatrix64F target = needsPost ? scratchK : gramian;

        fillScaledKernel(target, x1, x2, nR, nC, useWeights);

        // Optional post-processing (projection/normalization)
        if (needsPost) {
            postProcess(target, gramian, x1, nR); // pass x1 as default h-source
        }
    }

    private void fillScaledKernel(DenseMatrix64F out, double[] x1, double[] x2,
                                  int nR, int nC, boolean useWeights) {
        final double scale = kernel.getScale();
        for (int i = 0; i < nR; ++i) {
            final double xi = x1[i];                // TODO: generalize to multi-D row
            final double wi = useWeights ? w1[i] : 1.0;
            for (int j = 0; j < nC; ++j) {
                final double xj = x2[j];
                final double wj = useWeights ? w2[j] : 1.0;
                double v = kernel.getUnscaledCovariance(xi, xj) * wi * wj;
                out.add(i, j, scale * v );
            }
        }
    }

    protected void precomputeWeights(double[] x, int n, double[] out) {
        for (int i = 0; i < n; ++i) out[i] = weightFunction.getWeight(x[i]);
    }

    private void ensureCapacity(int nR, int nC) {
        if (scratchK == null || scratchK.numRows != nR || scratchK.numCols != nC)
            scratchK = new DenseMatrix64F(nR, nC);
        if (w1 == null || w1.length < nR) w1 = new double[nR];
        if (w2 == null || w2.length < nC) w2 = new double[nC];
    }

    // Subclass overrides only these two if needed.
    protected boolean needsPostProcess() { return false; }

    protected void postProcess(DenseMatrix64F inK, DenseMatrix64F out, double[] x1, int nR) {};

    public GaussianProcessKernel getKernel() { return kernel; }

    public DesignMatrix getDesignMatrix1() { return design1; }

    public DesignMatrix getDesignMatrix2() { return design2; }

    public WeightFunction getWeightFunction() {return weightFunction;}

    protected static DesignMatrix makeDesignMatrixFromWeights(RandomField.WeightProvider weights) {

        return new DesignMatrix("weights", false) {

            @Override
            public double getParameterValue(int row, int col) {
                throw new RuntimeException("Not yet implemented");
            }

            @Override
            public int getDimension() {
                return weights.getDimension();
            }

            @Override
            public int getRowDimension() {
                return weights.getDimension();
            }

            @Override
            public int getColumnDimension() {
                return 1;
            }

            @Override
            public Parameter getParameter(int column) {
                throw new IllegalArgumentException("Not allowed");
            }
        };
    }
}