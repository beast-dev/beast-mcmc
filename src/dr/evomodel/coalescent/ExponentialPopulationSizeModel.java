/*
 * ConstantPopulationModel.java
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

package dr.evomodel.coalescent;

import dr.evomodelxml.coalescent.demographicmodel.ExponentialGrowthModelParser;
import dr.inference.model.Parameter;

/**
 * Implements an exponential growth population size model.
 *
 * @author Andrew Rambaut
 */
public class ExponentialPopulationSizeModel extends PopulationSizeModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ExponentialPopulationSizeModel(Parameter N0Parameter, Parameter rateParameter, Type units) {

        this(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, N0Parameter, rateParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialPopulationSizeModel(String name, Parameter N0Parameter, Parameter rateParameter, Type units) {

        super(name, N0Parameter, units);

        this.rateParameter = rateParameter;
        addVariable(rateParameter);

        rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
    }

    public double getGrowthRate() {
        return rateParameter.getParameterValue(0);
    }

    /**
     * Return the population size function as an anonymous class...
     * @return the function
     */
    public PopulationSizeFunction getPopulationSizeFunction() {
        if (populationSizeFunction == null) {
            populationSizeFunction = new PopulationSizeFunction() {
                @Override
                public double getLogDemographic(double t) {
                    return getLogN0() + (-t * getGrowthRate());
                }

                @Override
                public double getIntegral(double startTime, double finishTime) {
                    double r = getGrowthRate();
                    if (r == 0.0) {
                        return (finishTime - startTime) / getN0();
                    } else {
                        return (Math.exp(finishTime * r) - Math.exp(startTime * r)) / getN0() / r;
                    }
                }

                @Override
                public Type getUnits() {
                    return ExponentialPopulationSizeModel.this.getUnits();
                }

                @Override
                public void setUnits(Type units) {
                    throw new UnsupportedOperationException();
                }
            };
        }

        return populationSizeFunction;
    }

    //
    // protected stuff
    //

    private final Parameter rateParameter;

    private PopulationSizeFunction populationSizeFunction = null;
}
