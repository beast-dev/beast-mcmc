/*
 * TransformedMultivariateHamiltonianMonteCarloOperatorParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.*;
import dr.util.Transform;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TransformedMultivariateHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {

    final static String TRANSFORMED_MULTIVARIATE_HMC = "transformedMultivariateHamiltonianMonteCarlo";

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions,
                                                    MassPreconditioner preconditioner, MassPreconditionScheduler.Type schedulerType,
                                                    ReversibleHMCProvider reversibleHMCprovider) {

        return new TransformedMultivariateHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner, schedulerType);
    }


    @Override
    public String getParserDescription() {
        return "Returns a Hamiltonian Monte Carlo transition kernel with dynamic mask on (only) transformed space";
    }

    @Override
    public Class getReturnType() {
        return TransformedMultivariateHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return TRANSFORMED_MULTIVARIATE_HMC;
    }

}
