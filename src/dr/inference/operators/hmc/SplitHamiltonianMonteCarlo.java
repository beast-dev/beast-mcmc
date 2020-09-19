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

    public SplitHamiltonianMonteCarlo(ReversibleHMCProvider inner,
                                      ReversibleHMCProvider outer, double stepSize,
                                      double relativeScale) {

        // TODO Call providers: inner and outer

        this.inner = inner; //todo: better names. Now A = HZZ, B = HMC
        this.outer = outer;
        this.dimA = inner.getInitialPosition().length;
        this.dimB = outer.getInitialPosition().length;

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
        outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        inner.reversiblePositionMomentumUpdate(positionA, momentumA, direction, relativeScale * time);
        outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        //3:merge the position and momentum, update position and momentum
        updateMergedVector(positionA, positionB, position);
        updateMergedVector(momentumA, momentumB, momentum);

        // TODO Avoid all these unncessary copies via

//        WrappedVector positionA = new WrappedVector.View(position, 0, dimA);
//        WrappedVector momentumA = new WrappedVector.View(momentum, 0, dimA);
//
//        WrappedVector positionB = new WrappedVector.View(position, dimA, dimB);
//        WrappedVector momentumB = new WrappedVector.View(momentum, dimA, dimB);
//
//        outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
//        inner.reversiblePositionMomentumUpdate(positionA, momentumA, direction, relativeScale * time);
//        outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        // done


    }


    @Override
    public double[] getInitialPosition() {

        double[] jointPosition = new double[dimA + dimB];
        System.arraycopy(inner.getInitialPosition(), 0, jointPosition, 0, dimA);
        System.arraycopy(outer.getInitialPosition(), 0, jointPosition, dimA, dimB);
        return jointPosition;
    }

    @Override
    public double getParameterLogJacobian() {
        return inner.getParameterLogJacobian() + outer.getParameterLogJacobian();
    }

    @Deprecated
    private WrappedVector mergeWrappedVector(WrappedVector vectorA, WrappedVector vectorB) {
        double[] buffer = new double[dimA + dimB];
        System.arraycopy(vectorA.getBuffer(), 0, buffer, 0, dimA);
        System.arraycopy(vectorB.getBuffer(), 0, buffer, dimA, dimB);
        return new WrappedVector.Raw(buffer);
    }

    @Deprecated
    private void splitWrappedVector(WrappedVector vectorAB, double[] bufferA, double[] bufferB) {
        System.arraycopy(vectorAB.getBuffer(), 0, bufferA, 0, dimA);
        System.arraycopy(vectorAB.getBuffer(), dimA, bufferB, 0, dimB);
    }

    @Deprecated
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

        inner.setParameter(bufferA);
        outer.setParameter(bufferB);

        // TODO Maybe could remove extra copy (if rate-limiting) by passing WrappedVector or ReadableVector
    }

    @Override
    public WrappedVector drawMomentum() {
        return mergeWrappedVector(inner.drawMomentum(), outer.drawMomentum());
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
//        if (reversibleHMCProviderB.getLogLikelihood() != reversibleHMCProviderA.getLogLikelihood()){
//            System.exit(-1);
//        }
        double[] momentumAbuffer = new double[dimA];
        double[] momentumBbuffer = new double[dimB];
        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);  // TODO Use WrappedVector.View instead
        //todo: better solution. Now only part B has a prior.
        return inner.getJointProbability(new WrappedVector.Raw(momentumAbuffer))
                + outer.getJointProbability(new WrappedVector.Raw(momentumBbuffer))
                - inner.getLogLikelihood();
    }

    @Override
    public double getLogLikelihood() {
        //todo: better solution. Now only part B has a prior. A only has likelihood.
        return inner.getLogLikelihood();
    }

    @Override
    public double getKineticEnergy(ReadableVector momentum) {

        double[] bufferA = new double[dimA];
        double[] bufferB = new double[dimB];

        System.arraycopy(((WrappedVector) momentum).getBuffer(), 0, bufferA, 0, dimA);
        System.arraycopy(((WrappedVector) momentum).getBuffer(), dimA, bufferB, 0, dimB);

        ReadableVector momentumA = new WrappedVector.Raw(bufferA);
        ReadableVector momentumB = new WrappedVector.Raw(bufferB);

        // TODO Use WrappedVector.View (or make a ReadableVector.View) instead

        return inner.getKineticEnergy(momentumA) + outer.getKineticEnergy(momentumB);
    }

    @Override
    public double getStepSize() { //todo:
        return 0.001;
    } //todo: tuning.

    private int dimA;
    private int dimB;
    private double stepSize;
    private double relativeScale;
    private ReversibleHMCProvider inner;
    private ReversibleHMCProvider outer;
}
