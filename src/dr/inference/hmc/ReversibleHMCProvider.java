package dr.inference.hmc;

import dr.math.matrixAlgebra.WrappedVector;

public interface ReversibleHMCProvider {
    void updatePositionAfterMap(WrappedVector position, WrappedVector momentum, double time);
    void updatePositionAfterMap(double[] position, WrappedVector momentum, double time);
}
