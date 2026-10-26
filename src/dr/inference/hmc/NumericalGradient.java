/*
 * NumericalGradient.java
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

import dr.inference.model.Bounds;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Andy Magee
 */
    public class NumericalGradient implements GradientWrtParameterProvider {

    public NumericalGradient(Likelihood likelihood, Parameter parameter) {
        this.parameter = parameter;
        this.likelihood = likelihood;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public GradientType getGradientType() {
        return GradientType.NUMERICAL;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    /**
     * Central differences with the step NumericalDerivative uses (the step size ratio times max(|x|, 1)), shortened to
     * stay inside the parameter's bounds: a positive rate close to 0, for example, would otherwise be evaluated below
     * its lower bound, where the density is NaN or -Inf. At a bound, a one-sided difference.
     */
    @Override
    public double[] getGradientLogDensity() {
        final double[] x = parameter.getParameterValues();
        final double[] gradient = new double[x.length];
        final Bounds<Double> bounds = boundsOf(parameter);
        final double ratio = getNumericGradientStepSize();

        for (int i = 0; i < x.length; ++i) {
            final double lower = bounds != null ? bounds.getLowerLimit(i) : Double.NEGATIVE_INFINITY;
            final double upper = bounds != null ? bounds.getUpperLimit(i) : Double.POSITIVE_INFINITY;
            final double room = Math.min(x[i] - lower, upper - x[i]);
            double h = ratio * Math.max(Math.abs(x[i]), 1.0);
            if (room > 0.0) {
                h = Math.min(h, room / 2.0);
                gradient[i] = (valueAt(i, x[i] + h) - valueAt(i, x[i] - h)) / (2.0 * h);
            } else if (x[i] - lower <= 0.0) {
                gradient[i] = (valueAt(i, x[i] + h) - valueAt(i, x[i])) / h;
            } else {
                gradient[i] = (valueAt(i, x[i]) - valueAt(i, x[i] - h)) / h;
            }
            parameter.setParameterValueQuietly(i, x[i]);
        }
        parameter.fireParameterChangedEvent();
        return gradient;
    }

    // null for a parameter without bounds, whose getBounds() throws
    private static Bounds<Double> boundsOf(Parameter parameter) {
        try {
            return parameter.getBounds();
        } catch (NullPointerException e) {
            return null;
        }
    }

    private double valueAt(int index, double value) {
        parameter.setParameterValueQuietly(index, value);
        parameter.fireParameterChangedEvent();
        return likelihood.getLogLikelihood();
    }

    private final Likelihood likelihood;
    private final Parameter parameter;
}
