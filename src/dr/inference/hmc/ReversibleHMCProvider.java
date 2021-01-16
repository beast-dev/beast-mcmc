package dr.inference.hmc;

import dr.inference.model.Parameter;
import dr.math.AdaptableCovariance;
import dr.math.matrixAlgebra.Lanczos;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
/**
 * @author Zhenyu Zhang
 */

public interface ReversibleHMCProvider {

    static double getMinEigValueLanczos(Parameter parameter, AdaptableCovariance sampleCov) {

        ReadableMatrix scmArray = sampleCov.getCovariance();
        double[] eigenvalues = Lanczos.eigen(scmArray, parameter.getDimension());

        if (eigenvalues.length < parameter.getDimension()) {
            throw new RuntimeException("called getMinEigValueSCM too early!");
        }

        System.err.println("largest eigenvalue is " + eigenvalues[0] + "smallest is " + eigenvalues[parameter.getDimension() - 1]);
        return eigenvalues[parameter.getDimension() - 1];
    }

    void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, WrappedVector gradient,
                                          int direction, double time);

    double[] getInitialPosition();

    double getParameterLogJacobian();

    int getNumGradientEvent();

    int getNumBoundaryEvent();

    Transform getTransform();

    GradientWrtParameterProvider getGradientProvider();

    void setParameter(double[] position);

    WrappedVector drawMomentum();

    double getJointProbability(WrappedVector momentum);

    double getLogLikelihood();

    double getKineticEnergy(ReadableVector momentum);

    double getStepSize();

    double getMinEigValueSCM();

    int getReversibleUpdateCount();

    boolean shouldUpdateSCM();
}
