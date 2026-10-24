/*
 * ExpansionModel.java
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

package dr.evomodel.coalescent.demographicmodel;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.Expansion;
import dr.evomodelxml.coalescent.demographicmodel.ExpansionModelParser;
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
 */
public class ExpansionModel extends DemographicModel implements Citable, DifferentiableDemographicModel {

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

    @Override
    public boolean isDifferentiableWrtParameter(Parameter parameter) {
        return parameter == N0Parameter || parameter == N1Parameter || parameter == growthRateParameter;
    }

    // arguments of Expansion: N0, r, N1, with N1 = p N0 where N1Parameter holds the proportion p
    @Override
    public double[] getArgumentDerivatives(Parameter parameter) {
        if (!isDifferentiableWrtParameter(parameter)) {
            throw new IllegalArgumentException("Not a parameter of " + getId() + ": " + parameter.getId());
        }
        final double N0 = N0Parameter.getParameterValue(0);
        final double proportion = N1Parameter.getParameterValue(0);

        double[] derivatives = new double[3];
        if (parameter == N0Parameter) {
            derivatives[0] += 1.0;
            derivatives[2] += proportion;
        }
        if (parameter == growthRateParameter) {
            derivatives[1] += getGrowthRateDerivative();
        }
        if (parameter == N1Parameter) {
            derivatives[2] += N0;
        }
        return derivatives;
    }

    // d r / d parameter, where the parameter is either the growth rate or the doubling time ln 2 / r
    private double getGrowthRateDerivative() {
        if (usingGrowthRate) {
            return 1.0;
        }
        final double doublingTime = growthRateParameter.getParameterValue(0);
        return -Math.log(2) / (doublingTime * doublingTime);
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
