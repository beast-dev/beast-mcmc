package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */

public class LogisticRegression extends GeneralizedLinearModel {


	public LogisticRegression(Parameter dependentParam) { //, Parameter independentParam, DesignMatrix designMatrix) {
		super(dependentParam);//, independentParam, designMatrix);
	}


	protected double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient) {
		return 0;  // todo
	}

	protected double calculateLogLikelihood(double[] beta) {
		// logLikelihood calculation for logistic regression
		throw new RuntimeException("Not yet implemented for optimization");
	}

	public boolean requiresScale() {
		return false;
	}

	protected double calculateLogLikelihood() {
		// logLikelihood calculation for logistic regression
		double logLikelihood = 0;

		double[] xBeta = getXBeta();

		for (int i = 0; i < N; i++) {
			// for each "pseudo"-datum
			logLikelihood += dependentParam.getParameterValue(i) * xBeta[i]
					- Math.log(1.0 + Math.exp(xBeta[i]));

		}
		return logLikelihood;
	}


	public boolean confirmIndependentParameters() {
		// todo -- check that independent parameters \in {0,1} only
		return true;
	}
}
