/*
 * ExpansionModel.java
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
import dr.evolution.coalescent.Expansion;
import dr.evomodelxml.coalescent.ExpansionModelParser;
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
 * @version $Id: ExpansionModel.java,v 1.5 2005/05/24 20:25:57 rambaut Exp $
 */
public class ExpansionModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ExpansionModel(Parameter N0Parameter, Parameter N1Parameter,
                          Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        this(ExpansionModelParser.EXPANSION_MODEL, N0Parameter, N1Parameter, growthRateParameter, units, usingGrowthRate);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExpansionModel(String name, Parameter N0Parameter, Parameter N1Parameter,
                          Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        super(name);

        expansion = new Expansion(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.N1Parameter = N1Parameter;
        addVariable(N1Parameter);
        N1Parameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.usingGrowthRate = usingGrowthRate;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        double N0 = N0Parameter.getParameterValue(0);
        double N1 = N1Parameter.getParameterValue(0);
        double growthRate = growthRateParameter.getParameterValue(0);

        if (usingGrowthRate) {
            expansion.setGrowthRate(growthRate);
        } else {
            double doublingTime = growthRate;
            growthRate = Math.log(2) / doublingTime;
            expansion.setDoublingTime(doublingTime);
        }

        expansion.setN0(N0);
        expansion.setProportion(N1);

        return expansion;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter N1Parameter = null;
    Parameter growthRateParameter = null;
    Expansion expansion = null;
    boolean usingGrowthRate = true;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Expansion Coalescent";
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
