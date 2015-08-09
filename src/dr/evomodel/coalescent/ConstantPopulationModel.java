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
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.model.Parameter;

/**
 * A wrapper for ConstantPopulation.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ConstantPopulationModel.java,v 1.10 2005/05/24 20:25:57 rambaut Exp $
 */
public class ConstantPopulationModel extends DemographicModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ConstantPopulationModel(Parameter N0Parameter, Type units) {

        this(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, N0Parameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstantPopulationModel(String name, Parameter N0Parameter, Type units) {

        super(name);

        constantPopulation = new ConstantPopulation(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        setUnits(units);
    }

    // general functions

    public DemographicFunction getDemographicFunction() {
        constantPopulation.setN0(N0Parameter.getParameterValue(0));
        return constantPopulation;
    }

    //
    // protected stuff
    //

    private Parameter N0Parameter;
    private ConstantPopulation constantPopulation = null;
}
