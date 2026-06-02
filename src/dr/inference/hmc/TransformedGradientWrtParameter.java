/*
 * TransformedGradientWrtParameter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ReciprocalLikelihood;
import dr.inference.model.TransformedParameter;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */

public class TransformedGradientWrtParameter implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradient;
    private final TransformedParameter parameter;
    private final ReciprocalLikelihood reciprocalLikelihood;
    private final boolean includeJacobian;
    private final boolean inverse;

    public TransformedGradientWrtParameter(GradientWrtParameterProvider gradient,
                                           TransformedParameter parameter,
                                           ReciprocalLikelihood reciprocalLikelihood,
                                           boolean includeJacobian,
                                           boolean inverse) {
        this.gradient = gradient;
        this.parameter = parameter;
        this.reciprocalLikelihood = reciprocalLikelihood;
        this.includeJacobian = includeJacobian;
        this.inverse = inverse;
    }
    @Override
    public Likelihood getLikelihood() {
        if (reciprocalLikelihood != null) {
            return reciprocalLikelihood;
        }
        return gradient.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return parameter.getUntransformedParameter();
    }

    @Override
    public int getDimension() {
        return parameter.getUntransformedParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] transformedGradient = gradient.getGradientLogDensity();
        double[] untransformedValues = parameter.getParameterUntransformedValues();

        Transform transform = parameter.getTransform();
        if (inverse) {
            transform = transform.inverseTransform();
        }

        double[] untransformedGradient;
        if (transform instanceof Transform.MultivariableTransform) {
            Transform.MultivariableTransform multivariableTransform = (Transform.MultivariableTransform) transform;
            untransformedGradient = multivariableTransform.updateGradientLogDensity(transformedGradient, untransformedValues,
                    0, untransformedValues.length);

            if (!includeJacobian) {
                throw new RuntimeException("Not yet implemented");
            }
        } else {

            double[] transformedValues = parameter.getParameterValues();
            untransformedGradient = new double[transformedGradient.length];
            for (int i = 0; i < untransformedGradient.length; ++i) {
                untransformedGradient[i] = transform.updateGradientLogDensity(transformedGradient[i], transformedValues[i]);
            }

            if (!includeJacobian) {
                for (int i = 0; i < untransformedGradient.length; ++i) {
                    untransformedGradient[i] -= transform.gradientLogJacobianInverse(parameter.getParameterValue(i));
                }
            }
        }

        if (reciprocalLikelihood != null) {
            for (int i = 0; i < untransformedGradient.length; ++i) {
                untransformedGradient[i] *= -1;
            }
        }

        return untransformedGradient;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, null);
    }
}
