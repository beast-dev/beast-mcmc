package dr.math.distributions.gp;

import dr.inference.model.DesignMatrix;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Filippo Monti
 * //
 * This class orthogonalizes the kernel against the linear subspace
 */
public class OrthogonalComplementProjectedBasisDimension extends BasisDimension {

    public OrthogonalComplementProjectedBasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2, WeightFunction weightFunction) {
        super(kernel, design1, design2, weightFunction);
    }
    public OrthogonalComplementProjectedBasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2) {
        this(kernel, design1, design2, null);
    }
    public OrthogonalComplementProjectedBasisDimension(GaussianProcessKernel kernel, DesignMatrix design) {
        this(kernel, design, design);
    }

    @Override
    protected boolean needsPostProcess() { return true; }

    @Override
    protected void postProcess(DenseMatrix64F inK, DenseMatrix64F out, double[] x1, int n) {
        final double[] v = new double[n]; // scratch
        // Adding the orthogonalized contribution of inK into 'out'
        computeKPerpAddInto(inK, out, x1, v);
//        DenseMatrix64F kPerp = computeKPerp(inK, x1); // slower approach
//        for(int i=0;i<n;i++) {
//            for (int j = 0; j < n; j++) {
//                out.add(i, j, kPerp.get(i, j));
//            }
//        }
    }

    /**
     * Adds K_perp (orthogonalized K) into 'out' without overwriting existing values.
     * Computes:
     *   out += K - (v h^T + h v^T)/s + ((h^T v)/s^2) h h^T,
     * with v = K h and s = h^T h.
     */
    private static void computeKPerpAddInto(DenseMatrix64F K, DenseMatrix64F out,
                                            double[] h, double[] v) {
        final int n = K.numRows;
        if (K.numCols != n || out.numRows != n || out.numCols != n)
            throw new IllegalArgumentException("K and out must be square n×n");

        // s = h^T h
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += h[i] * h[i];
        if (s <= 0.0) throw new IllegalArgumentException("h must be non-zero");
        final double invS = 1.0 / s;

        // v = K h
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            final int row = i * n;
            for (int j = 0; j < n; ++j) sum += K.data[row + j] * h[j];
            v[i] = sum;
        }
        // hv = h^T v
        double hv = 0.0;
        for (int i = 0; i < n; ++i) hv += h[i] * v[i];
        final double beta = hv * invS * invS;

        final double[] kd = K.data;
        final double[] od = out.data;

        // out += K - invS*(v h^T + h v^T) + beta*(h h^T)
        for (int i = 0; i < n; ++i) {
            final double vi = v[i];
            final double hi = h[i];
            final int row = i * n;
            for (int j = 0; j < n; ++j) {
                od[row + j] += kd[row + j]
                        - invS * (vi * h[j] + hi * v[j])
                        + beta * (hi * h[j]);
            }
        }

        // enforcing symmetry:
        // for (int i=0;i<n;i++) for (int j=i+1;j<n;j++){
        //   double a = 0.5*(out.get(i,j)+out.get(j,i));
        //   out.set(i,j,a); out.set(j,i,a);
        // }
    }


    private static DenseMatrix64F computeKPerp(DenseMatrix64F K, double[] h) {
        final int n = K.numRows;
        if (K.numCols != n) throw new IllegalArgumentException("K must be square");
        if (h.length != n)   throw new IllegalArgumentException("h length must match K dimension");

        // s = h^T h
        double s = 0.0;
        for (double hi : h) s += hi * hi;
        if (s <= 0.0) throw new IllegalArgumentException("h must be non-zero (h^T h > 0)");

        // v = K h
        final double[] v = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            int rowStart = i * n;
            for (int j = 0; j < n; j++) {
                sum += K.data[rowStart + j] * h[j];
            }
            v[i] = sum;
        }

        // hv = h^T v
        double hv = 0.0;
        for (int i = 0; i < n; i++) hv += h[i] * v[i];

        final double alpha = 1.0 / s;           // for (v h^T + h v^T)/s
        final double beta  = hv / (s * s);      // for (hv/s^2) h h^T

        // K_perp = K - alpha*(v h^T) - alpha*(h v^T) + beta*(h h^T)
        DenseMatrix64F Kperp = K.copy();
        double[] kd = Kperp.data;

        for (int i = 0; i < n; i++) {
            final double vi = v[i];
            final double hi = h[i];
            int rowStart = i * n;
            for (int j = 0; j < n; j++) {
                kd[rowStart + j] += -alpha * (vi * h[j] + hi * v[j]) + beta * (hi * h[j]);
            }
        }

        // TODO symmetrize to counter tiny round-off (since K should be symmetric)
        return Kperp;
    }
}