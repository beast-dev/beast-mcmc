/*
 * ExponentialConstantModel.java
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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExpConstant;
import dr.evomodelxml.coalescent.ExponentialConstantModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Exponential growth followed by constant size.
 *
 * @author Matthew Hall
 */
public class ExponentialConstantModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ExponentialConstantModel(Parameter N0Parameter,
                                    Parameter growthRateParameter,
                                    Parameter transitionTimeParameter,
                                    Type units) {

        this(ExponentialConstantModelParser.EXPONENTIAL_CONSTANT_MODEL,
                N0Parameter,
                growthRateParameter,
                transitionTimeParameter,
                units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialConstantModel(String name, Parameter N0Parameter,
                                    Parameter growthRateParameter,
                                    Parameter transitionTimeParameter,
                                    Type units) {

        super(name);

        exponentialConstant = new ExpConstant(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.transitionTimeParameter = transitionTimeParameter;
        addVariable(transitionTimeParameter);
        transitionTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        exponentialConstant.setN0(N0Parameter.getParameterValue(0));

        exponentialConstant.setGrowthRate(growthRateParameter.getParameterValue(0));

        exponentialConstant.setTransitionTime(transitionTimeParameter.getParameterValue(0));

        return exponentialConstant;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter growthRateParameter = null;
    Parameter transitionTimeParameter = null;
    ExpConstant exponentialConstant = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Exponential-Constant Coalescent";
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