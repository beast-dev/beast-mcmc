package dr.inference.hmc;

import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
/**
 * @author Zhenyu Zhang
 */

public interface ReversibleHMCProvider {

    void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, WrappedVector gradient,
                                          int direction, double time);
    void providerUpdatePreconditioning();

    double[] getInitialPosition();

    double getParameterLogJacobian();

    int getNumGradientEvent();// todo: make another interface to record these counts.

    int getNumBoundaryEvent();

    double[] getMask();

    Transform getTransform();

    GradientWrtParameterProvider getGradientProvider();

    void setParameter(double[] position);

    WrappedVector drawMomentum();

    double getJointProbability(WrappedVector momentum);

    double getLogLikelihood();

    double getKineticEnergy(ReadableVector momentum);

    double getStepSize();
}
