/*
 * RandomFieldGradient.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.RandomFieldDistribution;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Pratyusa Datta
 * @author Filippo Monti
 * @author Xiang Ji
 */

public class RandomFieldGradient implements GradientWrtParameterProvider, Reportable {

    private final RandomField randomField;
    private final Parameter parameter;
    private final RandomFieldDistribution distribution;
    private final GradientProvider provider;

    public RandomFieldGradient(RandomField randomField,
                               Parameter parameter) {

        this.randomField = randomField;
        this.parameter = parameter;
        this.distribution = randomField.getDistribution();

        if (parameter == randomField.getField()) {
            provider = new GradientProvider() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return distribution.getGradientLogDensity(x);
                }
            };
        } else {
            provider = distribution.getGradientWrt(parameter);
        }
    }
    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, null);
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

        // TODO Can cache here
        return provider.getGradientLogDensity(randomField.getField().getParameterValues());
    }
}