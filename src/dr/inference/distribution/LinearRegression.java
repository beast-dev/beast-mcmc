/*
 * LinearRegression.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
@Deprecated // GLM stuff is now in inference.glm - this is here for backwards compatibility temporarily
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
