package dr.inference.distribution;

import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 2, 2007
 * Time: 2:43:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class LogisticRegression extends GeneralizedLinearModel {


	public LogisticRegression(Parameter dependentParam, Parameter independentParam, DesignMatrix designMatrix) {
		super(dependentParam, independentParam, designMatrix);
	}


	protected double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient) {
		return 0;  //To change body of implemented methods use File | Settings | File Templates.
	}

	protected double calculateLogLikelihood(double[] beta) {
		// logLikelihood calculation for logistic regression
		double logLikelihood = 0;

		final int K = beta.length;
		final int N = dependentParam.getDimension();

		for (int i = 0; i < N; i++) {
			// for each "pseudo"-datum
			double xBeta = 0;
			for (int k = 0; k < K; k++) {
				xBeta += designMatrix.getParameterValue(i, k) * beta[k];
			}

			logLikelihood += dependentParam.getParameterValue(i) * xBeta
					- Math.log(1.0 + Math.exp(xBeta));

		}
		return logLikelihood;
	}

	protected boolean requiresScale() {
		return false;
	}

	protected double calculateLogLikelihood() {
		// logLikelihood calculation for logistic regression
		double logLikelihood = 0;

		final int K = independentParam.getDimension();
		final int N = dependentParam.getDimension();

		for (int i = 0; i < N; i++) {
			// for each "pseudo"-datum
			double xBeta = 0;
			for (int k = 0; k < K; k++) {
				xBeta += designMatrix.getParameterValue(i, k) * independentParam.getParameterValue(k);
			}

			logLikelihood += dependentParam.getParameterValue(i) * xBeta
					- Math.log(1.0 + Math.exp(xBeta));

		}
		return logLikelihood;
	}


	public boolean confirmIndependentParameters() {
		// todo -- check that independent parameters \in {0,1} only
		return true;
	}
}
