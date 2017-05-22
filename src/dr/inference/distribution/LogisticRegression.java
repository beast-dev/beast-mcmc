/*
 * LogisticRegression.java
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
 * @author Marc A. Suchard
 */
@Deprecated // GLM stuff is now in inference.glm - this is here for backwards compatibility temporarily
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
