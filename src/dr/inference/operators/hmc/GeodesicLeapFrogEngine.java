package dr.inference.operators.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;


public class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

    private final MatrixParameterInterface matrixParameter;

    GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                           MassPreconditioner preconditioning, double[] mask) {
        super(parameter, instabilityHandler, preconditioning, mask);
        this.matrixParameter = (MatrixParameterInterface) parameter;
    }

    @Override
    public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

        //TODO:
    }

    @Override
    public void updatePosition(double[] position, WrappedVector momentum,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

        //TODO:

    }

    @Override
    public void projectMomentum(WrappedVector momentum) {
        //TODO:
    }
}

