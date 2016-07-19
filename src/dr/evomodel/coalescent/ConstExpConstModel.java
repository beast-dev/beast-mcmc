/*
 * ConstExpConstModel.java
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

import dr.evolution.coalescent.ConstExpConst;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodelxml.coalescent.ConstExpConstModelParser;
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
public class ConstExpConstModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ConstExpConstModel(Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter,
                              Parameter epochParameter, boolean useNumericalIntegrator, Type units) {

        this(ConstExpConstModelParser.CONST_EXP_CONST_MODEL, N0Parameter, N1Parameter, growthRateParameter, timeParameter, epochParameter, useNumericalIntegrator, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstExpConstModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter,
                              Parameter epochParameter, boolean useNumericalIntegrator, Type units) {

        super(name);

        if (N1Parameter != null && growthRateParameter != null) {
            throw new RuntimeException("Only one of N1 and growthRate can be specified");
        }

        constExpConst = new ConstExpConst(
                (N1Parameter != null ? ConstExpConst.Parameterization.ANCESTRAL_POPULATION_SIZE : ConstExpConst.Parameterization.GROWTH_RATE),
                useNumericalIntegrator, units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.N1Parameter = N1Parameter;
        this.growthRateParameter = growthRateParameter;

        if (N1Parameter != null) {
            addVariable(N1Parameter);
            N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        } else {
            addVariable(growthRateParameter);
            growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        }

        this.timeParameter = timeParameter;
        addVariable(timeParameter);
        timeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.epochParameter = epochParameter;
        addVariable(epochParameter);
        epochParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        constExpConst.setEpochTime(epochParameter.getParameterValue(0));
        constExpConst.setN0(N0Parameter.getParameterValue(0));
        constExpConst.setTime1(timeParameter.getParameterValue(0));

        if (N1Parameter != null) {
            constExpConst.setN1(N1Parameter.getParameterValue(0));
        } else {
            constExpConst.setGrowthRate(growthRateParameter.getParameterValue(0));
        }

        return constExpConst;
    }

    //
    // protected stuff
    //

    private final Parameter N0Parameter;
    private final Parameter N1Parameter;
    private final Parameter growthRateParameter;
    private final Parameter timeParameter;
    private final Parameter epochParameter;
    private final ConstExpConst constExpConst;


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Constant-Exponential-Constant Coalescent";
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