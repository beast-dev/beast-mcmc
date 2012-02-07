package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class LinearRegression extends GeneralizedLinearModel {

	private static final double normalizingConstant = -0.5 * Math.log(2 * Math.PI);

          private boolean logTransform = false;

        public double[] getTransformedDependentParameter() {
            double[] y = dependentParam.getParameterValues();
            if (logTransform) {
                for(int i=0; i<y.length; i++)
                    y[i] = Math.log(y[i]);
            }
            return y;
        }

	protected double calculateLogLikelihood() {
		double logLikelihood = 0;
		double[] xBeta = getXBeta();
		double[] precision = getScale();
                    double[] y = getTransformedDependentParameter();
              
		for (int i = 0; i < N; i++) {    // assumes that all observations are independent given fixed and random effects
                              if (logTransform)
                                  logLikelihood -= y[i]; // Jacobian
			logLikelihood += 0.5 * Math.log(precision[i]) - 0.5 * (y[i] - xBeta[i]) * (y[i] - xBeta[i]) * precision[i];

		}
		return N * normalizingConstant + logLikelihood;
	}

	public LinearRegression(Parameter dependentParam, boolean logTransform) { //, Parameter independentParam, DesignMatrix designMatrix) {
		super(dependentParam); //, independentParam, designMatrix);
		System.out.println("Constructing a linear regression model");
                    this.logTransform = logTransform;
	}

	protected double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient) {
                    throw new RuntimeException("Optimization not yet implemented.");
	}

	public boolean requiresScale() {
		return true;
	}

	protected double calculateLogLikelihood(double[] beta) {
		throw new RuntimeException("Optimization not yet implemented.");
	}

	protected boolean confirmIndependentParameters() {
		return true;
	}
}
