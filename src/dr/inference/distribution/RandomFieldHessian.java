/*
 * RandomFieldHessian.java
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

package dr.inference.distribution;

import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.RandomFieldDistribution;

public class RandomFieldHessian implements HessianWrtParameterProvider {

    private final RandomField randomField;
    private final Parameter parameter;
    private final RandomFieldDistribution distribution;

    public RandomFieldHessian(RandomField randomField,
                              Parameter parameter) {
        this.randomField = randomField;
        this.parameter = parameter;
        this.distribution = randomField.getDistribution();

        if (parameter != randomField.getField()) {
            throw new IllegalArgumentException("Random-field Hessian is only implemented for the field parameter");
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return randomField;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return distribution.getGradientLogDensity(randomField.getField().getParameterValues());
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return distribution.getDiagonalHessianLogDensity(randomField.getField().getParameterValues());
    }

    @Override
    public double[][] getHessianLogDensity() {
        return distribution.getHessianLogDensity(randomField.getField().getParameterValues());
    }
}
