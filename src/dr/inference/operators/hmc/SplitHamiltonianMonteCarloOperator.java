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
    private final int updateRelativeScaleDelay = 1000;
    private final int updateRelativeScaleFrequency = 1000;
    private ReversibleHMCProvider inner;
    private ReversibleHMCProvider outer;
    protected final Parameter parameter;

    int dimInner;
    int dimOuter;

    private int nSteps;
    private int nSplitOuter;


    private int gradientCheckCount;
    private double gradientCheckTolerance;


    public SplitHamiltonianMonteCarloOperator(double weight, ReversibleHMCProvider inner,
                                              ReversibleHMCProvider outer, Parameter parameter,
                                              double stepSize,
                                              double relativeScale, int nSteps, int nSplitOuter,
                                              int gradientCheckCount, double gradientCheckTolerance) {

        setWeight(weight);
        this.inner = inner;
        this.outer = outer;
        dimInner = inner.getInitialPosition().length;
        dimOuter = outer.getInitialPosition().length;
        this.parameter = parameter;

        this.stepSize = stepSize;
        this.relativeScale = relativeScale;
        this.nSteps = nSteps;
        this.nSplitOuter = nSplitOuter;

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
        if (parameter.getDimension() != dimInner + dimOuter) {
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

        if (shouldUpdateRelativeScale()){
            relativeScale = Math.sqrt(inner.getMinEigValueSCM()) / Math.sqrt(outer.getMinEigValueSCM());
        }

        double[] positionInnerbuffer = inner.getInitialPosition();
        double[] positionOuterbuffer = outer.getInitialPosition();

        WrappedVector positionInner = new WrappedVector.Raw(positionInnerbuffer);
        WrappedVector positionOuter = new WrappedVector.Raw(positionOuterbuffer);

        WrappedVector momentumInner = inner.drawMomentum();
        WrappedVector momentumOuter = outer.drawMomentum();

        WrappedVector gradientInner = new WrappedVector.Raw(inner.getGradientProvider().getGradientLogDensity());
        WrappedVector gradientOuter = new WrappedVector.Raw(outer.getGradientProvider().getGradientLogDensity());

        final double prop =
                inner.getKineticEnergy(momentumInner) + outer.getKineticEnergy(momentumOuter) + inner.getParameterLogJacobian() + outer.getParameterLogJacobian();

        for (int i = 0; i < nSteps; i++) {
            for (int j = 0; j < nSplitOuter; j++) {
                outer.reversiblePositionMomentumUpdate(positionOuter, momentumOuter, gradientOuter,1, .5 * stepSize / nSplitOuter);
            }
            inner.reversiblePositionMomentumUpdate(positionInner, momentumInner, gradientInner, 1,
                    relativeScale * stepSize);
            updateOuterGradient(gradientOuter);
            for (int j = 0; j < nSplitOuter; j++) {
                outer.reversiblePositionMomentumUpdate(positionOuter, momentumOuter, gradientOuter,1, .5 * stepSize / nSplitOuter);
            }
        }

        final double res =
                inner.getKineticEnergy(momentumInner) + outer.getKineticEnergy(momentumOuter) + inner.getParameterLogJacobian() + outer.getParameterLogJacobian();

        return prop - res;
    }

    public void updateOuterGradient(WrappedVector gradient) {
        double[] buffer = outer.getGradientProvider().getGradientLogDensity();
        for (int i = 0; i < buffer.length; i++) {
            gradient.set(i, buffer[i]);
        }
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
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, WrappedVector gradient, int direction,
                                                 double time) {
        if (shouldUpdateRelativeScale()){
            relativeScale = Math.sqrt(inner.getMinEigValueSCM()) / Math.sqrt(outer.getMinEigValueSCM());
        }

        double[] positionInnerbuffer = new double[dimInner];
        double[] positionOuterbuffer = new double[dimOuter];

        double[] momentumAbuffer = new double[dimInner];
        double[] momentumBbuffer = new double[dimOuter];

        //1:split the position
        splitWrappedVector(position, positionInnerbuffer, positionOuterbuffer);

        //2:split the momentum
        splitWrappedVector(momentum, momentumAbuffer, momentumBbuffer);

        WrappedVector positionInner = new WrappedVector.Raw(positionInnerbuffer);
        WrappedVector positionOuter = new WrappedVector.Raw(positionOuterbuffer);
        WrappedVector momentumA = new WrappedVector.Raw(momentumAbuffer);
        WrappedVector momentumB = new WrappedVector.Raw(momentumBbuffer);

        //2:update them
        for (int i = 0; i < nSplitOuter; i++) {
            outer.reversiblePositionMomentumUpdate(positionOuter, momentumB, gradient, direction, .5 * time / nSplitOuter);
        }
        inner.reversiblePositionMomentumUpdate(positionInner, momentumA, gradient, direction, relativeScale * time);
        updateOuterGradient(gradient);
        for (int i = 0; i < nSplitOuter; i++) {
            outer.reversiblePositionMomentumUpdate(positionOuter, momentumB, gradient, direction, .5 * time / nSplitOuter);
        }
        //3:merge the position and momentum, update position and momentum
        updateMergedVector(positionInner, positionOuter, position);
        updateMergedVector(momentumA, momentumB, momentum);
    }

    public double[] jointTransformInverse(double[] argument) {
        double[] jointUntransformedPosition = new double[dimInner + dimOuter];
        double[] transformedPositionB = new double[dimOuter];

        System.arraycopy(argument, dimInner, transformedPositionB, 0, dimOuter);
        double[] unTransformedPositionB = outer.getTransform().inverse(transformedPositionB, 0, dimOuter);

        System.arraycopy(argument, 0, jointUntransformedPosition, 0, dimInner);
        System.arraycopy(unTransformedPositionB, 0, jointUntransformedPosition, dimInner, dimOuter);

        return jointUntransformedPosition;
    }

    public double transformGetLogJacobian(double[] untransformedValue) {
        double[] untransformedPositionB = new double[dimOuter];
        System.arraycopy(untransformedValue, dimInner, untransformedPositionB, 0, dimOuter);
        return outer.getTransform().getLogJacobian(untransformedPositionB, 0, dimOuter);
    }

    public double[] transformupdateGradientLogDensity(double[] analyticalGradientOriginal, Parameter parameter) {

        double[] analyticalGradientOriginalInner = new double[dimInner];
        double[] analyticalGradientOriginalOuter = new double[dimOuter];
        double[] parameterValueOuter = new double[dimOuter];
        double[] updatedGradientJoint = new double[dimInner + dimOuter];
        System.arraycopy(analyticalGradientOriginal, 0, analyticalGradientOriginalInner, 0, dimInner);
        System.arraycopy(analyticalGradientOriginal, dimInner, analyticalGradientOriginalOuter, 0, dimOuter);
        System.arraycopy(parameter.getParameterValues(), dimInner, parameterValueOuter, 0, dimOuter);

        double[] updatedGradientOuter =
                outer.getTransform().updateGradientLogDensity(analyticalGradientOriginalOuter,
                        parameterValueOuter, 0, dimOuter);
        System.arraycopy(analyticalGradientOriginalInner, 0, updatedGradientJoint, 0, dimInner);
        System.arraycopy(updatedGradientOuter, 0, updatedGradientJoint, dimInner, dimOuter);
        return updatedGradientJoint;
    }

    @Override
    public double[] getInitialPosition() {

        double[] jointPosition = new double[dimInner + dimOuter];
        System.arraycopy(inner.getInitialPosition(), 0, jointPosition, 0, dimInner);
        System.arraycopy(outer.getInitialPosition(), 0, jointPosition, dimInner, dimOuter);
        return jointPosition;
    }

    @Override
    public double getParameterLogJacobian() {
        return inner.getParameterLogJacobian() + outer.getParameterLogJacobian();
    }

    @Override
    public int getNumGradientEvent() {
        return 0;
    }

    @Override
    public int getNumBoundaryEvent() {
        return 0;
    }

    @Override
    public Transform getTransform() {
        return null;
    }

    @Override
    public GradientWrtParameterProvider getGradientProvider() {
        return outer.getGradientProvider();
    }

    private double[] mergeGradient() {
        double[] jointPosition = new double[dimInner + dimOuter];
        System.arraycopy(inner.getGradientProvider().getGradientLogDensity(), 0, jointPosition, 0,
                dimInner);
        System.arraycopy(outer.getGradientProvider().getGradientLogDensity(), 0, jointPosition, dimInner
                , dimOuter);
        return jointPosition;
    }


    @Override
    public void setParameter(double[] position) {
        double[] bufferA = new double[dimInner];
        double[] bufferB = new double[dimOuter];

        System.arraycopy(position, 0, bufferA, 0, dimInner);
        System.arraycopy(position, dimInner, bufferB, 0, dimOuter);

        inner.setParameter(bufferA);
        outer.setParameter(bufferB);
    }

    @Override
    public WrappedVector drawMomentum() {
        return mergeWrappedVector(inner.drawMomentum(), outer.drawMomentum());
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
        double[] momentumAbuffer = new double[dimInner];
        double[] momentumBbuffer = new double[dimOuter];
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
        double[] bufferA = new double[dimInner];
        double[] bufferB = new double[dimOuter];

        System.arraycopy(((WrappedVector) momentum).getBuffer(), 0, bufferA, 0, dimInner);
        System.arraycopy(((WrappedVector) momentum).getBuffer(), dimInner, bufferB, 0, dimOuter);

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

    @Override
    public double getMinEigValueSCM() {
        return 1;//todo
    }

    private boolean shouldUpdateRelativeScale(){
        //use updateRelativeScaleDelay updateRelativeScaleFrequency getCount
        return true;
    }

    private WrappedVector mergeWrappedVector(WrappedVector lhs, WrappedVector rhs) {
        double[] buffer = new double[lhs.getDim() + rhs.getDim()];
        System.arraycopy(lhs.getBuffer(), 0, buffer, 0, lhs.getDim());
        System.arraycopy(rhs.getBuffer(), 0, buffer, lhs.getDim(), rhs.getDim());
        return new WrappedVector.Raw(buffer);
    }

    private void splitWrappedVector(WrappedVector vectorJoint, double[] bufferLHS, double[] bufferRHS) {
        System.arraycopy(vectorJoint.getBuffer(), 0, bufferLHS, 0, dimInner);
        System.arraycopy(vectorJoint.getBuffer(), dimInner, bufferRHS, 0, dimOuter);
    }

    private void updateMergedVector(WrappedVector vectorLHS, WrappedVector vectorRHS, WrappedVector vectorJoint) {
        for (int i = 0; i < dimInner + dimOuter; i++) {
            if (i < dimInner) {
                vectorJoint.set(i, vectorLHS.get(i));
            } else {
                vectorJoint.set(i, vectorRHS.get(i - dimInner));
            }
        }
    }
}
