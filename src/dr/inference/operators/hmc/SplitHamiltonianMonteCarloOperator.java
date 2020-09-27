package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.GeneralOperator;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

/**
 * @author Zhenyu Zhang
 * @author Aki Nishimura
 * @author Marc A. Suchard
 */

public class SplitHamiltonianMonteCarloOperator extends AbstractAdaptableOperator implements GeneralOperator,
        ReversibleHMCProvider {

    private double stepSize;
    private double relativeScale;
    private ReversibleHMCProvider inner;
    private ReversibleHMCProvider outer;
    protected final Parameter parameter;

    int dimA;
    int dimB;

    private int nSteps;
    private int innerSteps;


    private int gradientCheckCount;
    private double gradientCheckTolerance;


    public SplitHamiltonianMonteCarloOperator(double weight, ReversibleHMCProvider inner,
                                              ReversibleHMCProvider outer, Parameter parameter,
                                              double stepSize,
                                              double relativeScale, int nSteps, int innerSteps,
                                              int gradientCheckCount, double gradientCheckTolerance) {

        setWeight(weight);
        this.inner = inner; //todo: better names. Now A = HZZ, B = HMC
        this.outer = outer;
        dimA = inner.getInitialPosition().length;
        dimB = outer.getInitialPosition().length;
        this.parameter = parameter;

        this.stepSize = stepSize;
        this.relativeScale = relativeScale;
        this.nSteps = nSteps;
        this.innerSteps = innerSteps;

        this.gradientCheckCount = gradientCheckCount;
        this.gradientCheckTolerance = gradientCheckTolerance;
    }

    public double doOperation(Likelihood joint) {
        if (getCount() < gradientCheckCount) {
            checkGradient(joint);
        }
        return mergedUpdate();
    }

    void checkGradient(final Likelihood joint) {
        if (parameter.getDimension() != dimA + dimB) {
            throw new RuntimeException("Unequal dimensions");
        }

        MultivariateFunction numeric = new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                if (!anyTransform()) {

                    ReadableVector.Utils.setParameter(argument, parameter);
                    return joint.getLogLikelihood();
                } else {

                    double[] untransformedValue = jointTransformInverse(argument);
                    ReadableVector.Utils.setParameter(untransformedValue, parameter);
                    return joint.getLogLikelihood() - transformGetLogJacobian(untransformedValue);
                }
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return parameter.getBounds().getLowerLimit(n);
            }

            @Override
            public double getUpperBound(int n) {
                return parameter.getBounds().getUpperLimit(n);
            }
        };

        double[] analyticalGradientOriginal = mergeGradient();

        double[] restoredParameterValue = parameter.getParameterValues();

        if (!anyTransform()) {

            double[] numericGradientOriginal = NumericalDerivative.gradient(numeric, parameter.getParameterValues());

            if (!MathUtils.isClose(analyticalGradientOriginal, numericGradientOriginal, gradientCheckTolerance)) {

                String sb = "Gradients do not match:\n" +
                        "\tAnalytic: " + new WrappedVector.Raw(analyticalGradientOriginal) + "\n" +
                        "\tNumeric : " + new WrappedVector.Raw(numericGradientOriginal) + "\n";
                throw new RuntimeException(sb);
            }

        } else {

            double[] transformedParameter = getInitialPosition();
            double[] numericGradientTransformed = NumericalDerivative.gradient(numeric, transformedParameter);

            double[] analyticalGradientTransformed = transformupdateGradientLogDensity(analyticalGradientOriginal,
                    parameter);

            if (!MathUtils.isClose(analyticalGradientTransformed, numericGradientTransformed, gradientCheckTolerance)) { //todo: read in the tolerance.
                String sb = "Transformed Gradients do not match:\n" +
                        "\tAnalytic: " + new WrappedVector.Raw(analyticalGradientTransformed) + "\n" +
                        "\tNumeric : " + new WrappedVector.Raw(numericGradientTransformed) + "\n" +
                        "\tParameter : " + new WrappedVector.Raw(parameter.getParameterValues()) + "\n" +
                        "\tTransformed Parameter : " + new WrappedVector.Raw(transformedParameter) + "\n";
                throw new RuntimeException(sb);
            }
        }

        ReadableVector.Utils.setParameter(restoredParameterValue, parameter);
    }

    private boolean anyTransform() {
        return (inner.getParameterLogJacobian() == 0 && outer.getParameterLogJacobian() == 0) ? false : true;
    }

    private double mergedUpdate() {

        double[] positionAbuffer = inner.getInitialPosition();
        double[] positionBbuffer = outer.getInitialPosition();

        WrappedVector positionA = new WrappedVector.Raw(positionAbuffer);
        WrappedVector positionB = new WrappedVector.Raw(positionBbuffer);

        WrappedVector momentumA = inner.drawMomentum();
        WrappedVector momentumB = outer.drawMomentum();

        final double prop =
                inner.getKineticEnergy(momentumA) + outer.getKineticEnergy(momentumB) + inner.getParameterLogJacobian() + outer.getParameterLogJacobian();

        for (int i = 0; i < nSteps; i++) {
            for (int j = 0; j < innerSteps; j++) {
                outer.reversiblePositionMomentumUpdate(positionB, momentumB, 1, .5 * stepSize / innerSteps);
            }
            inner.reversiblePositionMomentumUpdate(positionA, momentumA, 1,
                    relativeScale * stepSize);
            for (int j = 0; j < innerSteps; j++) {
                outer.reversiblePositionMomentumUpdate(positionB, momentumB, 1, .5 * stepSize / innerSteps);
            }
        }

        final double res =
                inner.getKineticEnergy(momentumA) + outer.getKineticEnergy(momentumB) + inner.getParameterLogJacobian() + outer.getParameterLogJacobian();

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

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, int direction,
                                                 double time) {

        double[] positionAbuffer = new double[dimA];
        double[] positionBbuffer = new double[dimB];

        double[] momentumAbuffer = new double[dimA];
        double[] momentumBbuffer = new double[dimB];

        //1:split the position
        splitWrappedVector(position, positionAbuffer, positionBbuffer);

        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);

        WrappedVector positionA = new WrappedVector.Raw(positionAbuffer);
        WrappedVector positionB = new WrappedVector.Raw(positionBbuffer);
        WrappedVector momentumA = new WrappedVector.Raw(momentumAbuffer);
        WrappedVector momentumB = new WrappedVector.Raw(momentumBbuffer);

        //2:update them
        for (int i = 0; i < innerSteps; i++) {
            outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, .5 * time / innerSteps);
        }
        inner.reversiblePositionMomentumUpdate(positionA, momentumA, direction, relativeScale * time);
        for (int i = 0; i < innerSteps; i++) {
            outer.reversiblePositionMomentumUpdate(positionB, momentumB, direction, .5 * time / innerSteps);
        }
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

    }

    public double[] jointTransformInverse(double[] argument) {
        double[] jointUntransformedPosition = new double[dimA + dimB];
        double[] transformedPositionB = new double[dimB];

        System.arraycopy(argument, dimA, transformedPositionB, 0, dimB);
        double[] unTransformedPositionB = outer.getTransform().inverse(transformedPositionB, 0, dimB);

        System.arraycopy(argument, 0, jointUntransformedPosition, 0, dimA);
        System.arraycopy(unTransformedPositionB, 0, jointUntransformedPosition, dimA, dimB);

        return jointUntransformedPosition;
    }

    public double transformGetLogJacobian(double[] untransformedValue) {
        double[] untransformedPositionB = new double[dimB];
        System.arraycopy(untransformedValue, dimA, untransformedPositionB, 0, dimB);
        return outer.getTransform().getLogJacobian(untransformedPositionB, 0, dimB);
    }

    public double[] transformupdateGradientLogDensity(double[] analyticalGradientOriginal, Parameter parameter) {

        double[] analyticalGradientOriginalA = new double[dimA];
        double[] analyticalGradientOriginalB = new double[dimB];
        double[] parameterValueB = new double[dimB];
        double[] updatedGradientAB = new double[dimA + dimB];
        System.arraycopy(analyticalGradientOriginal, 0, analyticalGradientOriginalA, 0, dimA);
        System.arraycopy(analyticalGradientOriginal, dimA, analyticalGradientOriginalB, 0, dimB);
        System.arraycopy(parameter.getParameterValues(), dimA, parameterValueB, 0, dimB);

        double[] updatedGradientB =
                outer.getTransform().updateGradientLogDensity(analyticalGradientOriginalB,
                        parameterValueB, 0, dimB);
        System.arraycopy(analyticalGradientOriginalA, 0, updatedGradientAB, 0, dimA);
        System.arraycopy(updatedGradientB, 0, updatedGradientAB, dimA, dimB);
        return updatedGradientAB;
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

    @Override
    public Transform getTransform() {
        return null;
    }

    @Override
    public GradientWrtParameterProvider getGradientProvider() {
        return null;
    }

    private double[] mergeGradient() {
        double[] jointPosition = new double[dimA + dimB];
        System.arraycopy(inner.getGradientProvider().getGradientLogDensity(), 0, jointPosition, 0,
                dimA);
        System.arraycopy(outer.getGradientProvider().getGradientLogDensity(), 0, jointPosition, dimA
                , dimB);
        return jointPosition;
    }


    @Override
    public void setParameter(double[] position) {
        double[] bufferA = new double[dimA];
        double[] bufferB = new double[dimB];

        System.arraycopy(position, 0, bufferA, 0, dimA);
        System.arraycopy(position, dimA, bufferB, 0, dimB);

        inner.setParameter(bufferA);
        outer.setParameter(bufferB);
    }

    @Override
    public WrappedVector drawMomentum() {
        return mergeWrappedVector(inner.drawMomentum(), outer.drawMomentum());
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
        double[] momentumAbuffer = new double[dimA];
        double[] momentumBbuffer = new double[dimB];
        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);
        //todo: better solution. Now only part B has a prior.
        return inner.getJointProbability(new WrappedVector.Raw(momentumAbuffer)) + outer.getJointProbability(new WrappedVector.Raw(momentumBbuffer)) - inner.getLogLikelihood();
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
        // TODO Use ReadableVector.View instead
//        return inner.getKineticEnergy(new ReadableVector.View(momentum, 0, dimA)) +
//                outer.getKineticEnergy(new ReadableVector.View(momentum, dimA, dimB));

        return inner.getKineticEnergy(momentumA) + outer.getKineticEnergy(momentumB);
    }

    @Override
    public double getStepSize() { //todo:
        return stepSize;
    } //todo: tuning.

    private WrappedVector mergeWrappedVector(WrappedVector lhs, WrappedVector rhs) {
        double[] buffer = new double[lhs.getDim() + rhs.getDim()];
        System.arraycopy(lhs.getBuffer(), 0, buffer, 0, lhs.getDim());
        System.arraycopy(rhs.getBuffer(), 0, buffer, lhs.getDim(), rhs.getDim());
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
}
