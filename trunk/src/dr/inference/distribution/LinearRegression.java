package dr.inference.distribution;

import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 2, 2007
 * Time: 3:58:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class LinearRegression extends GeneralizedLinearModel {

    protected double calculateLogLikelihood() {
        return 0;
    }

    public LinearRegression(Parameter dependentParam, Parameter independentParam, DesignMatrix designMatrix) {
        super(dependentParam, independentParam, designMatrix);
    }

    protected double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient) {
        return 0;
    }


    protected boolean requiresScale() {
        return true;
    }

    protected double calculateLogLikelihood(double[] beta) {
        double logLikelihood = 0;

        final int K = beta.length;
        final int N = dependentParam.getDimension();

        for (int i = 0; i < N; i++) {
            // for each "pseudo"-datum
            double xBeta = 0;
            for (int k = 0; k < K; k++) {
                xBeta += designMatrix.getParameterValue(i, k) * beta[k];
            }

//            logLikelihood += dependentParam.getParameterValue(i) * xBeta
//                    - Math.log(1.0 + Math.exp(xBeta));
//                                               todo
        }
        return logLikelihood;
    }

    protected boolean confirmIndependentParameters() {
        return true;
    }
}
