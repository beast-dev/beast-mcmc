/*
 * TransformedMultivariateHamiltonianMonteCarloOperator.java
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


package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TransformedMultivariateHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    private final Parameter maskParameter;

    public TransformedMultivariateHamiltonianMonteCarloOperator(AdaptationMode mode,
                                                                double weight,
                                                                GradientWrtParameterProvider gradientProvider,
                                                                Parameter parameter,
                                                                Transform transform,
                                                                Parameter maskParameter,
                                                                Options runtimeOptions,
                                                                MassPreconditioner preconditioner,
                                                                MassPreconditionScheduler.Type preconditionSchedulerType) {
        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner, preconditionSchedulerType);

        this.maskParameter = maskParameter;

    }

    @Override
    protected LeapFrogEngine constructLeapFrogEngine(Transform transform) {
        return new MaskedMultivariateTransform(parameter, getDefaultInstabilityHandler(), preconditioning, maskParameter, (Transform.MultivariateTransform) transform);
    }

    @Override
    public String getReport() {
        return null;
    }

    class MaskedMultivariateTransform extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

        private final Transform.MultivariateTransform transform;
        private final Parameter mask;


        MaskedMultivariateTransform(Parameter parameter,
                                    InstabilityHandler instabilityHandler,
                                    MassPreconditioner preconditioning,
                                    Parameter mask,
                                    Transform.MultivariateTransform transform) {
            super(parameter, instabilityHandler, preconditioning, mask.getParameterValues());
            this.transform = transform;
            this.mask = mask;
            setMaskUntransformedSpace();
        }

        private void setMaskUntransformedSpace() {
            for (int i = 0; i < super.mask.length; i++) {
                super.mask[i] = 1.0;
            }
        }
    }

}
