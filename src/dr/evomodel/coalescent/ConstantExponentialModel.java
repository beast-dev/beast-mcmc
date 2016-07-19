/*
 * ConstantExponentialModel.java
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

import dr.evolution.coalescent.ConstExponential;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodelxml.coalescent.ConstantExponentialModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Exponential growth from a constant ancestral population size.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ConstantExponentialModel.java,v 1.8 2005/10/28 02:49:17 alexei Exp $
 */
public class ConstantExponentialModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ConstantExponentialModel(Parameter N0Parameter, Parameter timeParameter,
                                    Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        this(ConstantExponentialModelParser.CONSTANT_EXPONENTIAL_MODEL, N0Parameter, timeParameter, growthRateParameter, units, usingGrowthRate);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstantExponentialModel(String name, Parameter N0Parameter, Parameter timeParameter,
                                    Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        super(name);

        constExponential = new ConstExponential(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.timeParameter = timeParameter;
        addVariable(timeParameter);
        timeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.usingGrowthRate = usingGrowthRate;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        double time = timeParameter.getParameterValue(0);
        double N0 = N0Parameter.getParameterValue(0);
        double growthRate = growthRateParameter.getParameterValue(0);

        if (!usingGrowthRate) {
            double doublingTime = growthRate;
            growthRate = Math.log(2) / doublingTime;
        }

        constExponential.setGrowthRate(growthRate);
        constExponential.setN0(N0);
        constExponential.setN1(N0 * Math.exp(-time * growthRate));

        return constExponential;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter timeParameter = null;
    Parameter growthRateParameter = null;
    ConstExponential constExponential = null;
    boolean usingGrowthRate = true;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Constant-Exponential Coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("OG", "Pybus"),
                    new Author("A", "Rambaut")
            },
            "GENIE: estimating demographic history from molecular phylogenies",
            2001,
            "Bioinformatics",
            18, 1404, 1405
    );
}
