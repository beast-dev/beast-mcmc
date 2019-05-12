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

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.wrightfisher.Population;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.model.Parameter;

/**
 * Implements an constant size population size model.
 *
 * @author Andrew Rambaut
 */
public class ConstantPopulationSizeModel extends PopulationSizeModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ConstantPopulationSizeModel(Parameter N0Parameter, boolean inLogSpace, Type units) {

        this(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, N0Parameter, inLogSpace, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstantPopulationSizeModel(String name, Parameter N0Parameter, boolean inLogSpace, Type units) {

        super(name, N0Parameter, units);
    }

    /**
     * Return the population size function as an anonymous class...
     * @return the function
     */
    public PopulationSizeFunction getPopulationSizeFunction() {
        return new PopulationSizeFunction() {
            @Override
            public double getLogDemographic(double t) {
                return getLogN0();
            }

            @Override
            public double getIntegral(double startTime, double finishTime) {
                return (finishTime - startTime) / getN0();
            }

            @Override
            public Type getUnits() {
                return ConstantPopulationSizeModel.this.getUnits();
            }

            @Override
            public void setUnits(Type units) {
                throw new UnsupportedOperationException();
            }
        };
    }

}
