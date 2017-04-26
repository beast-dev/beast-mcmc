package dr.inference.operators;

import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

/**
 * @author Max Tolkoff
 */
public class HMCOperator extends AbstractCoercableOperator{
    GradientProvider derivative;
    Parameter parameter;
    double stepSize;
    int nSteps;
    NormalDistribution drawDistribution;


    public HMCOperator(CoercionMode mode, double weight, GradientProvider derivative, Parameter parameter, double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);
        this.derivative = derivative;
        this.parameter = parameter;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));

    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "HMC Operator";
    }

    @Override
    public double doOperation() { //throws OperatorFailedException {
        double functionalStepSize = stepSize;

        double[] HMCDerivative = derivative.getGradientLogDensity();
        double[] momentum = new double[HMCDerivative.length];
        for (int i = 0; i < momentum.length; i++) {
            momentum[i] = (Double) drawDistribution.nextRandom();
        }

        double prop=0;
        for (int i = 0; i < momentum.length ; i++) {
            prop += momentum[i] * momentum[i] / (2 * Math.pow(drawDistribution.getSD(), 2));
        }


        for (int i = 0; i < momentum.length; i++) {
                momentum[i] = momentum[i] - functionalStepSize / 2 * HMCDerivative[i];
        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j < momentum.length; j++) {
                    parameter.setParameterValue(j, parameter.getParameterValue(j) + functionalStepSize * momentum[j] / (Math.pow(drawDistribution.getSD() , 2)));
            }
//            parameter.fireParameterChangedEvent();

            HMCDerivative = derivative.getGradientLogDensity();

            if(i != nSteps){

                for (int j = 0; j < momentum.length; j++) {
                    momentum[j] = momentum[j] - functionalStepSize / 2 * HMCDerivative[j];
                }
            }
        }

        for (int i = 0; i < momentum.length; i++) {
                momentum[i] = momentum[i] - functionalStepSize / 2 * HMCDerivative[i];
        }

        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res += momentum[i] * momentum[i] / (2 * Math.pow(drawDistribution.getSD(), 2));
        }
        return prop - res;
    }

    @Override
    public double getCoercableParameter() {
        return Math.log(stepSize);
    }

    @Override
    public void setCoercableParameter(double value) {
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }
}
