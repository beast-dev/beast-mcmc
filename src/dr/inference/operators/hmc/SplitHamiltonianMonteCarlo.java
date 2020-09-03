package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Zhenyu Zhang
 * @author Aki Nishimura
 * @author Marc A. Suchard
 */

public class SplitHamiltonianMonteCarlo implements ReversibleHMCProvider {

    public SplitHamiltonianMonteCarlo(ReversibleHMCProvider reversibleHMCProviderA,
                                      ReversibleHMCProvider reversibleHMCProviderB, double stepSize,
                                      double relativeScale) {

        this.reversibleHMCProviderA = reversibleHMCProviderA; //todo: better names. Now A = HZZ, B = HMC
        this.reversibleHMCProviderB = reversibleHMCProviderB;
        dimA = reversibleHMCProviderA.getInitialPosition().length;
        dimB = reversibleHMCProviderB.getInitialPosition().length;

        this.stepSize = stepSize;
        this.relativeScale = relativeScale;
    }

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, int direction,
                                                 double time) {

        double[] positionAbuffer = new double[dimA];
        double[] positionBbuffer = new double[dimB];

        double[] momentumAbuffer = new double[dimA];
        double[] momentumBbuffer = new double[dimB];

        //1:split the position (the order is (A,B) for position and momentum)
        splitWrappedVector(position, positionAbuffer, positionBbuffer);

        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);

        WrappedVector positionA = new WrappedVector.Raw(positionAbuffer);
        WrappedVector positionB = new WrappedVector.Raw(positionBbuffer);
        WrappedVector momentumA = new WrappedVector.Raw(momentumAbuffer);
        WrappedVector momentumB = new WrappedVector.Raw(momentumBbuffer);
        // WrappedVector (check offset)?

        //2:update them
        reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        reversibleHMCProviderA.reversiblePositionMomentumUpdate(positionA, momentumA, direction, relativeScale * time);
        reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        //3:merge the position and momentum, update position and momentum
        updateMergedVector(positionA, positionB, position);
        updateMergedVector(momentumA, momentumB, momentum);
    }


    @Override
    public double[] getInitialPosition() {

        double[] jointPosition = new double[dimA + dimB];
        System.arraycopy(reversibleHMCProviderA.getInitialPosition(), 0, jointPosition, 0, dimA);
        System.arraycopy(reversibleHMCProviderB.getInitialPosition(), 0, jointPosition, dimA, dimB);
        return jointPosition;
    }

    @Override
    public double getParameterLogJacobian() {
        return reversibleHMCProviderA.getParameterLogJacobian() + reversibleHMCProviderB.getParameterLogJacobian();
    }

    private WrappedVector mergeWrappedVector(WrappedVector vectorA, WrappedVector vectorB) {
        double[] buffer = new double[dimA + dimB];
        System.arraycopy(vectorA.getBuffer(), 0, buffer, 0, dimA);
        System.arraycopy(vectorB.getBuffer(), 0, buffer, dimA, dimB);
        return new WrappedVector.Raw(buffer);
    }

    private void splitWrappedVector(WrappedVector vectorAB, double[] bufferA, double[] bufferB) {
        System.arraycopy(vectorAB.getBuffer(), 0, bufferA, 0, dimA);
        System.arraycopy(vectorAB.getBuffer(), dimA, bufferB, 0, dimB);
    }

    private void updateMergedVector(WrappedVector vectorA, WrappedVector vectorB, WrappedVector vectorAB) {
        for (int i = 0; i < dimA + dimB; i++) {
            if (i < dimA) {
                vectorAB.set(i, vectorA.get(i));
            } else {
                vectorAB.set(i, vectorB.get(i - dimA));
            }
        }
    }

    @Override
    public void setParameter(double[] position) {
        double[] bufferA = new double[dimA];
        double[] bufferB = new double[dimB];

        System.arraycopy(position, 0, bufferA, 0, dimA);
        System.arraycopy(position, dimA, bufferB, 0, dimB);

        reversibleHMCProviderA.setParameter(bufferA);
        reversibleHMCProviderB.setParameter(bufferB);
    }

    @Override
    public WrappedVector drawMomentum() {
        return mergeWrappedVector(reversibleHMCProviderA.drawMomentum(), reversibleHMCProviderB.drawMomentum());
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
//        if (reversibleHMCProviderB.getLogLikelihood() != reversibleHMCProviderA.getLogLikelihood()){
//            System.exit(-1);
//        }
        double[] momentumAbuffer = new double[dimA];
        double[] momentumBbuffer = new double[dimB];
        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);
        //todo: better solution. Now only part B has a prior.
        return reversibleHMCProviderA.getJointProbability(new WrappedVector.Raw(momentumAbuffer)) + reversibleHMCProviderB.getJointProbability(new WrappedVector.Raw(momentumBbuffer)) - reversibleHMCProviderA.getLogLikelihood();
    }

    @Override
    public double getLogLikelihood() {
        //todo: better solution. Now only part B has a prior. A only has likelihood.
        return reversibleHMCProviderA.getLogLikelihood();
    }

    @Override
    public double getKineticEnergy(ReadableVector momentum) {

        double[] bufferA = new double[dimA];
        double[] bufferB = new double[dimB];

        System.arraycopy(((WrappedVector) momentum).getBuffer(), 0, bufferA, 0, dimA);
        System.arraycopy(((WrappedVector) momentum).getBuffer(), dimA, bufferB, 0, dimB);

        ReadableVector momentumA = new WrappedVector.Raw(bufferA);
        ReadableVector momentumB = new WrappedVector.Raw(bufferB);

        return reversibleHMCProviderA.getKineticEnergy(momentumA) + reversibleHMCProviderB.getKineticEnergy(momentumB);
    }

    @Override
    public double getStepSize() { //todo:
        return 0.001;
    } //todo: tuning.

    private int dimA;
    private int dimB;
    private double stepSize;
    private double relativeScale;
    private ReversibleHMCProvider reversibleHMCProviderA;
    private ReversibleHMCProvider reversibleHMCProviderB;
}
