package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.math.AdaptableCovariance;
import dr.math.matrixAlgebra.Lanczos;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * An interface that provide a travel time multiplier for split HMC. The multiplier is a function of the sample
 * covariances of
 * the inner and outer operators.
 *
 * @author Zhenyu Zhang
 */
public interface SplitHMCtravelTimeMultiplier {

    double getMultiplier();

    void updateSCM(AdaptableCovariance scm, double[] value, long iter);

    boolean shouldUpdateSCM(long iter);

    boolean shouldGetMultiplier(long iter);

    AdaptableCovariance getInnerCov();

    AdaptableCovariance getOuterCov();

    static SplitHMCtravelTimeMultiplier create(ReversibleHMCProvider inner, ReversibleHMCProvider outer,
                                               RSoptions rsOptions) {
        return new LeastEigenvalueRatioSCM(inner, outer, rsOptions);
    }

    class RSoptions {

        final int updateRSdelay;
        final int updateRSfrequency;
        final int getRSdelay;

        public RSoptions(int updateRSdelay, int updateRSfrequency, int getRSdelay) {
            this.updateRSdelay = updateRSdelay;
            this.updateRSfrequency = updateRSfrequency;
            this.getRSdelay = getRSdelay;
        }
    }
}

class LeastEigenvalueRatioSCM implements SplitHMCtravelTimeMultiplier {

    ReversibleHMCProvider inner;
    ReversibleHMCProvider outer;
    AdaptableCovariance adaptableCovarianceInner;
    AdaptableCovariance adaptableCovarianceOuter;
    double[] maskVectorInner;
    double[] maskVectorOuter;

    RSoptions rsOptions;

    public LeastEigenvalueRatioSCM(ReversibleHMCProvider inner, ReversibleHMCProvider outer, RSoptions rsOptions) {

        this.inner = inner;
        this.outer = outer;
        this.rsOptions = rsOptions;

        this.adaptableCovarianceInner = new AdaptableCovariance(inner.getInitialPosition().length);
        this.adaptableCovarianceOuter = new AdaptableCovariance(outer.getInitialPosition().length);

        this.maskVectorInner = inner.getMask();
        this.maskVectorOuter = outer.getMask();
    }

    @Override
    public double getMultiplier() {
        double minEigenInner = getMinEigValueLanczos(adaptableCovarianceInner, maskVectorInner);
        double minEigenOuter = getMinEigValueLanczos(adaptableCovarianceOuter, maskVectorOuter);
        return Math.sqrt(minEigenInner) / Math.sqrt(minEigenOuter);
    }

    static ReadableMatrix subsetByMask(ReadableMatrix adaptableCov, double[] mask) {

        int sumDim = 0;

        for (int i = 0; i < mask.length; i++) {
            sumDim += mask[i];
        }

        if (sumDim < 2) {
            throw new RuntimeException("not a matrix!");
        }

        double[][] subsetMat = new double[sumDim][sumDim];

        int m = 0;
        int n;
        boolean allZeroRow;

        for (int i = 0; i < adaptableCov.getMajorDim(); i++) {
            n = 0;
            allZeroRow = true;
            for (int j = 0; j < adaptableCov.getMajorDim(); j++) {
                if (mask[i] == 1 && mask[j] == 1) {
                    allZeroRow = false;
                    subsetMat[m][n] = adaptableCov.get(i, j);
                    n++;
                }
            }
            if (!allZeroRow) m++;
        }
        return new WrappedMatrix.ArrayOfArray(subsetMat);
    }

    static double getMinEigValueLanczos(AdaptableCovariance sampleCov, double[] mask) {

        ReadableMatrix scmArray = mask != null ? subsetByMask(sampleCov.getCovariance(), mask) :
                sampleCov.getCovariance();
        int dim = scmArray.getMajorDim();
        double[] eigenvalues = Lanczos.eigen(scmArray, dim);

        if (eigenvalues.length < dim) {
            throw new RuntimeException("called getMinEigValueSCM too early!");
        }

        //System.err.println("largest eigenvalue is " + eigenvalues[0] + "smallest is " + eigenvalues[dim - 1]);
        return eigenvalues[dim - 1];
    }

    @Override
    public void updateSCM(AdaptableCovariance scm, double[] parameterValue, long iter) {

        if (shouldUpdateSCM(iter)) {
            scm.update(new WrappedVector.Raw(parameterValue));
        }
    }

    @Override
    public AdaptableCovariance getInnerCov() {
        return adaptableCovarianceInner;
    }

    @Override
    public AdaptableCovariance getOuterCov() {
        return adaptableCovarianceOuter;
    }

    public boolean shouldUpdateSCM(long iter) {
        return ((rsOptions.updateRSfrequency > 0)
                && ((iter % rsOptions.updateRSfrequency == 0)
                && (iter > rsOptions.updateRSdelay)));
    }

    public boolean shouldGetMultiplier(long iter) {
        return (iter > rsOptions.getRSdelay);
    }
}
