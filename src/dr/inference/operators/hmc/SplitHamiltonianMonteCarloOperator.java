package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Likelihood;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.GeneralOperator;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

import java.util.Arrays;
/**
 * @author Zhenyu Zhang
 * @author Aki Nishimura
 * @author Marc A. Suchard
 */

public class SplitHamiltonianMonteCarloOperator extends AbstractAdaptableOperator implements GeneralOperator{//todo: make a ReversibleHMCProvider for NUTS use.

    private double stepSize;
    private double relativeScale;
    private ReversibleHMCProvider reversibleHMCProviderA;
    private ReversibleHMCProvider reversibleHMCProviderB;
    private int nSteps;
    private int innerSteps;

    public SplitHamiltonianMonteCarloOperator(double weight, ReversibleHMCProvider reversibleHMCProviderA,
                                              ReversibleHMCProvider reversibleHMCProviderB, double stepSize,
                                              double relativeScale, int nSteps, int innerSteps) {
        setWeight(weight);
        this.reversibleHMCProviderA = reversibleHMCProviderA; //todo: better names. Now A = HZZ, B = HMC
        this.reversibleHMCProviderB = reversibleHMCProviderB;

        this.stepSize = stepSize;
        this.relativeScale = relativeScale;
        this.nSteps = nSteps;
        this.innerSteps = innerSteps;
    }

    public double doOperation(Likelihood joint) {
        return mergedUpdate();
    }

    private double mergedUpdate() {

        double[] positionAbuffer = reversibleHMCProviderA.getInitialPosition();
        double[] positionBbuffer = reversibleHMCProviderB.getInitialPosition();

        WrappedVector positionA = new WrappedVector.Raw(positionAbuffer);
        WrappedVector positionB = new WrappedVector.Raw(positionBbuffer);

        WrappedVector momentumA = reversibleHMCProviderA.drawMomentum();
        WrappedVector momentumB = reversibleHMCProviderB.drawMomentum();

        final double prop =
                reversibleHMCProviderA.getKineticEnergy(momentumA) + reversibleHMCProviderB.getKineticEnergy(momentumB) + reversibleHMCProviderA.getParameterLogJacobian() + reversibleHMCProviderB.getParameterLogJacobian();

        for (int i = 0; i < nSteps; i++) {
            for (int j = 0; j < innerSteps; j++) {
                reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, 1, stepSize);
            }
            reversibleHMCProviderA.reversiblePositionMomentumUpdate(positionA, momentumA, 1,
                    relativeScale * stepSize);
            for (int j = 0; j < innerSteps; j++) {
                reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, 1, stepSize);
            }
        }

        final double res =
                reversibleHMCProviderA.getKineticEnergy(momentumA) + reversibleHMCProviderB.getKineticEnergy(momentumB) + reversibleHMCProviderA.getParameterLogJacobian() + reversibleHMCProviderB.getParameterLogJacobian();

        return prop - res;
    }


    @Override
    public double doOperation() {
        throw new RuntimeException("Should not be executed");
    }

    @Override
    public String getOperatorName() {
        return "Split HMC operator";
    }


    @Override
    protected void setAdaptableParameterValue(double value) {

    }

    @Override
    protected double getAdaptableParameterValue() {
        return 1;
    }

    @Override
    public double getRawParameter() {
        return 1;
    }

    @Override
    public String getAdaptableParameterName() {
        return null;
    }
}
