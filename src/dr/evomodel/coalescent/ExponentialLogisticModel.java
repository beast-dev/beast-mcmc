/*
 * ExponentialLogisticModel.java
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
import dr.evolution.coalescent.ExponentialLogistic;
import dr.evomodelxml.coalescent.ExponentialLogisticModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Exponential growth followed by Logistic growth.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class ExponentialLogisticModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ExponentialLogisticModel(Parameter N0Parameter,
                                    Parameter logisticGrowthParameter,
                                    Parameter logisticShapeParameter,
                                    Parameter exponentialGrowthParameter,
                                    Parameter transitionTimeParameter,
                                    double alpha, Type units) {

        this(ExponentialLogisticModelParser.EXPONENTIAL_LOGISTIC_MODEL,
                N0Parameter,
                logisticGrowthParameter,
                logisticShapeParameter,
                exponentialGrowthParameter,
                transitionTimeParameter,
                alpha, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialLogisticModel(String name, Parameter N0Parameter,
                                    Parameter logisticGrowthParameter,
                                    Parameter logisticShapeParameter,
                                    Parameter exponentialGrowthParameter,
                                    Parameter transistionTimeParameter,
                                    double alpha, Type units) {

        super(name);

        exponentialLogistic = new ExponentialLogistic(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.logisticGrowthParameter = logisticGrowthParameter;
        addVariable(logisticGrowthParameter);
        logisticGrowthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.logisticShapeParameter = logisticShapeParameter;
        addVariable(logisticShapeParameter);
        logisticShapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.exponentialGrowthParameter = exponentialGrowthParameter;
        addVariable(exponentialGrowthParameter);
        exponentialGrowthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.transistionTimeParameter = transistionTimeParameter;
        addVariable(transistionTimeParameter);
        transistionTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.alpha = alpha;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        exponentialLogistic.setN0(N0Parameter.getParameterValue(0));

        double r = logisticGrowthParameter.getParameterValue(0);
        exponentialLogistic.setGrowthRate(r);

        double r1 = exponentialGrowthParameter.getParameterValue(0);
        exponentialLogistic.setR1(r1);

        double t = transistionTimeParameter.getParameterValue(0);
        exponentialLogistic.setTime(t);

        // logisticGrowth.setShape(Math.exp(shapeParameter.getParameterValue(0)));
        exponentialLogistic.setTime50(logisticShapeParameter.getParameterValue(0));
        //exponentialLogistic.setShapeFromTimeAtAlpha(logisticShapeParameter.getParameterValue(0), alpha);

        return exponentialLogistic;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter logisticGrowthParameter = null;
    Parameter logisticShapeParameter = null;
    Parameter exponentialGrowthParameter = null;
    Parameter transistionTimeParameter = null;
    double alpha = 0.5;
    ExponentialLogistic exponentialLogistic = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Exponential-Logistic Coalescent";
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