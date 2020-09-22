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
        ReversibleHMCProvider {//todo: implement ReversibleHMCProvider.

    private double stepSize;
    private double relativeScale;
    private ReversibleHMCProvider reversibleHMCProviderA;
    private ReversibleHMCProvider reversibleHMCProviderB;
    protected final Parameter parameter;

    int dimA;
    int dimB;

    private int nSteps;
    private int innerSteps;


    private int gradientCheckCount;
    private double gradientCheckTolerance;


    public SplitHamiltonianMonteCarloOperator(double weight, ReversibleHMCProvider reversibleHMCProviderA,
                                              ReversibleHMCProvider reversibleHMCProviderB, Parameter parameter,
                                              double stepSize,
                                              double relativeScale, int nSteps, int innerSteps,
                                              int gradientCheckCount, double gradientCheckTolerance) {

        setWeight(weight);
        this.reversibleHMCProviderA = reversibleHMCProviderA; //todo: better names. Now A = HZZ, B = HMC
        this.reversibleHMCProviderB = reversibleHMCProviderB;
        dimA = reversibleHMCProviderA.getInitialPosition().length;
        dimB = reversibleHMCProviderB.getInitialPosition().length;
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
        return (reversibleHMCProviderA.getParameterLogJacobian() == 0 && reversibleHMCProviderB.getParameterLogJacobian() == 0) ? false : true;
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
            reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        }
        reversibleHMCProviderA.reversiblePositionMomentumUpdate(positionA, momentumA, direction, relativeScale * time);
        for (int i = 0; i < innerSteps; i++) {
            reversibleHMCProviderB.reversiblePositionMomentumUpdate(positionB, momentumB, direction, time);
        }
        //3:merge the position and momentum, update position and momentum
        updateMergedVector(positionA, positionB, position);
        updateMergedVector(momentumA, momentumB, momentum);
    }

    public double[] jointTransformInverse(double[] argument) {
        double[] jointUntransformedPosition = new double[dimA + dimB];
        double[] transformedPositionB = new double[dimB];

        System.arraycopy(argument, dimA, transformedPositionB, 0, dimB);
        double[] unTransformedPositionB = reversibleHMCProviderB.getTransform().inverse(transformedPositionB, 0, dimB);

        System.arraycopy(argument, 0, jointUntransformedPosition, 0, dimA);
        System.arraycopy(unTransformedPositionB, 0, jointUntransformedPosition, dimA, dimB);

        return jointUntransformedPosition;
    }

    public double transformGetLogJacobian(double[] untransformedValue) {
        double[] untransformedPositionB = new double[dimB];
        System.arraycopy(untransformedValue, dimA, untransformedPositionB, 0, dimB);
        return reversibleHMCProviderB.getTransform().getLogJacobian(untransformedPositionB, 0, dimB);
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
                reversibleHMCProviderB.getTransform().updateGradientLogDensity(analyticalGradientOriginalB,
                        parameterValueB, 0, dimB);
        System.arraycopy(analyticalGradientOriginalA, 0, updatedGradientAB, 0, dimA);
        System.arraycopy(updatedGradientB, 0, updatedGradientAB, dimA, dimB);
        return updatedGradientAB;
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
        System.arraycopy(reversibleHMCProviderA.getGradientProvider().getGradientLogDensity(), 0, jointPosition, 0,
                dimA);
        System.arraycopy(reversibleHMCProviderB.getGradientProvider().getGradientLogDensity(), 0, jointPosition, dimA
                , dimB);
        return jointPosition;
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
        return stepSize;
    } //todo: tuning.

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
}
