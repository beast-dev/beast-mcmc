package dr.inference.distribution;

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;

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

//              System.err.println("N = "+N);
//              System.exit(-1);

//              System.err.println("XBeta = "+new Vector(xBeta));
//              System.exit(-1);
                    double[] y = getTransformedDependentParameter();
              
		for (int i = 0; i < N; i++) {    // assumes that all observations are independent given fixed and random effects

//			double y = dependentParam.getParameterValue(i);
                              if (logTransform) {
//                                  y = Math.log(y);
                                  logLikelihood -= y[i]; // Jacobian
                              }
//                        System.err.println("y = "+y);
//                        System.err.println("p_Y = "+precision[i]);
			logLikelihood += 0.5 * Math.log(precision[i]) - 0.5 * (y[i] - xBeta[i]) * (y[i] - xBeta[i]) * precision[i];

		}
//              System.err.println("");
//              System.exit(-1);
		return N * normalizingConstant + logLikelihood;
	}

	public LinearRegression(Parameter dependentParam, boolean logTransform) { //, Parameter independentParam, DesignMatrix designMatrix) {
		super(dependentParam); //, independentParam, designMatrix);
		System.out.println("Constructing a linear regression model");
                    this.logTransform = logTransform;
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
