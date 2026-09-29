package dr.math.distributions.gp;

import dr.inference.model.DesignMatrix;
import org.ejml.data.DenseMatrix64F;
/**
 * @author Filippo Monti
 * subclass of BasisDimension that fixes the variance of the unscaled linear kernel to one
 */

public class UnitaryVarianceLinearBasisDimension extends BasisDimension {

    public UnitaryVarianceLinearBasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2, WeightFunction weightFunction) {
        super(kernel, design1, design2, weightFunction);
    }
    public UnitaryVarianceLinearBasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2) {
        this(kernel, design1, design2, null);
    }
    public UnitaryVarianceLinearBasisDimension(GaussianProcessKernel kernel, DesignMatrix design) {
        this(kernel, design, design);
    }

    @Override
    protected boolean needsPostProcess() { return true; }

    @Override
    protected void postProcess(DenseMatrix64F inK, DenseMatrix64F out, double[] x1, int n) {

        double normalizer = 0;
        for (int i = 0; i < n; i++) normalizer += x1[i] * x1[i];
        normalizer /= n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out.add(i, j, inK.get(i, j) / normalizer);
            }
        }
    }
}