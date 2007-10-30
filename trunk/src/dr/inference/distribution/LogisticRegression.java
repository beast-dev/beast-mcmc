package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 2, 2007
 * Time: 2:43:23 PM
 * To change this template use File | Settings | File Templates.
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

	protected boolean requiresScale() {
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
