package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class LinearRegression extends GeneralizedLinearModel {

	private static final double normalizingConstant = -0.5 * Math.log(2 * Math.PI);

	protected double calculateLogLikelihood() {
		double logLikelihood = 0;
		double[] xBeta = getXBeta();
		double[] precision = getScale();

		for (int i = 0; i < N; i++) {    // assumes that all observations are independent given fixed and random effects

			double y = dependentParam.getParameterValue(i);
			logLikelihood += 0.5 * Math.log(precision[i]) - 0.5 * (y - xBeta[i]) * (y - xBeta[i]) * precision[i];

		}
		return N * normalizingConstant + logLikelihood;
	}

	public LinearRegression(Parameter dependentParam) { //, Parameter independentParam, DesignMatrix designMatrix) {
		super(dependentParam); //, independentParam, designMatrix);
		System.out.println("Constructing a linear regression model");
	}

	protected double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient) {
		return 0;
	}


	protected boolean requiresScale() {
		return true;
	}

	protected double calculateLogLikelihood(double[] beta) {
		throw new RuntimeException("Optimization not yet implemented.");
	}

	protected boolean confirmIndependentParameters() {
		return true;
	}
}
