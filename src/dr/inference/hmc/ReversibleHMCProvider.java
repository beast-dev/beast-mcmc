package dr.inference.hmc;

import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
/**
 * @author Zhenyu Zhang
 */

public interface ReversibleHMCProvider {

    void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, int direction, double time);

    double[] getInitialPosition();

    double getParameterLogJacobian();

    void setParameter(double[] position);

    WrappedVector drawMomentum();

    double getJointProbability(WrappedVector momentum);

    double getKineticEnergy(ReadableVector momentum);

    double getStepSize();

}
