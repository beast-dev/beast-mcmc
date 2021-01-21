package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.math.AdaptableCovariance;
import dr.math.matrixAlgebra.Lanczos;
import dr.math.matrixAlgebra.ReadableMatrix;
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

    RSoptions rsOptions;

    public LeastEigenvalueRatioSCM(ReversibleHMCProvider inner, ReversibleHMCProvider outer, RSoptions rsOptions) {

        this.inner = inner;
        this.outer = outer;
        this.rsOptions = rsOptions;

        this.adaptableCovarianceInner = new AdaptableCovariance(inner.getInitialPosition().length);
        this.adaptableCovarianceOuter = new AdaptableCovariance(outer.getInitialPosition().length);
    }

    @Override
    public double getMultiplier() {
        double minEigenInner = getMinEigValueLanczos(inner.getGradientProvider().getParameter(),
                adaptableCovarianceInner);
        double minEigenOuter = getMinEigValueLanczos(outer.getGradientProvider().getParameter(),
                adaptableCovarianceInner);
        return Math.sqrt(minEigenInner) / Math.sqrt(minEigenOuter);
    }

    static double getMinEigValueLanczos(Parameter parameter, AdaptableCovariance sampleCov) {

        ReadableMatrix scmArray = sampleCov.getCovariance();
        double[] eigenvalues = Lanczos.eigen(scmArray, parameter.getDimension());

        if (eigenvalues.length < parameter.getDimension()) {
            throw new RuntimeException("called getMinEigValueSCM too early!");
        }

        System.err.println("largest eigenvalue is " + eigenvalues[0] + "smallest is " + eigenvalues[parameter.getDimension() - 1]);
        return eigenvalues[parameter.getDimension() - 1];
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
