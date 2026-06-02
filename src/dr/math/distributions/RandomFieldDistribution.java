/*
 * RandomFieldDistribution.java
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

package dr.math.distributions;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.*;

/**
 * @author Pratyusa Datta
 * @author Filippo Monti
 * @author Xiang Ji
 */
public abstract class RandomFieldDistribution extends AbstractModel implements
        ParametricMultivariateDistributionModel, GradientProvider, HessianProvider {

    public RandomFieldDistribution(String name) {
        super(name);
    }

    public abstract GradientProvider getGradientWrt(Parameter parameter);

    public double getIncrement(int i, Parameter field) {
        double[] mean = getMean();
        return (field.getParameterValue(i) - mean[i]) - (field.getParameterValue(i + 1) - mean[i + 1]);
    }
}
