/*
 * ConstantLogisticModel.java
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

import dr.evolution.coalescent.ConstLogistic;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodelxml.coalescent.ConstantLogisticModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Logistic growth from a constant ancestral population size.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ConstantLogisticModel.java,v 1.7 2005/04/11 11:24:39 alexei Exp $
 */
public class ConstantLogisticModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ConstantLogisticModel(Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, Type units) {

        this(ConstantLogisticModelParser.CONSTANT_LOGISTIC_MODEL, N0Parameter, N1Parameter, growthRateParameter, shapeParameter, alpha, units);
    }

    /**
     * Construct demographic model with default settings
     */
    private ConstantLogisticModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, Type units) {

        super(name);

        constLogistic = new ConstLogistic(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.N1Parameter = N1Parameter;
        addVariable(N1Parameter);
        N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.alpha = alpha;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        constLogistic.setN0(N0Parameter.getParameterValue(0));
        constLogistic.setN1(N1Parameter.getParameterValue(0));

        double r = growthRateParameter.getParameterValue(0);
        constLogistic.setGrowthRate(r);

        // AER 24/02/03
        // logisticGrowth.setShape(Math.exp(shapeParameter.getParameterValue(0)));

        // New parameterization of logistic shape to be the time at which the
        // population reached some proportion alpha:
        double C = ((1.0 - alpha) * Math.exp(-r * shapeParameter.getParameterValue(0))) / alpha;
        constLogistic.setShape(C);

        return constLogistic;
    }

    //
    // protected stuff
    //

    private Parameter N0Parameter = null;
    private Parameter N1Parameter = null;
    private Parameter growthRateParameter = null;
    private Parameter shapeParameter = null;
    private double alpha = 0.5;
    private ConstLogistic constLogistic = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Constant-Logistic Coalescent";
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
