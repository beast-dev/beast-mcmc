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
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Implements an exponential growth population size model.
 *
 * @author Andrew Rambaut
 */
public class LogisticPopulationSizeModel extends PopulationSizeModel implements Citable {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public LogisticPopulationSizeModel(Parameter logN0Parameter, Parameter rateParameter, Parameter shapeParameter, Type units) {

        this(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, logN0Parameter, rateParameter, shapeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public LogisticPopulationSizeModel(String name, Parameter logN0Parameter, Parameter rateParameter, Parameter shapeParameter, Type units) {

        super(name, logN0Parameter, units);

        this.rateParameter = rateParameter;
        addVariable(rateParameter);

        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);

        rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getGrowthRate() {
        return rateParameter.getParameterValue(0);
    }

    public double getShape() {
        return shapeParameter.getParameterValue(0);
    }

    /**
     * Return the population size function as an anonymous class...
     * @return the function
     */
    public PopulationSizeFunction getPopulationSizeFunction() {
        return new PopulationSizeFunction() {
            @Override
            public double getLogDemographic(double t) {
                throw new UnsupportedOperationException("not implemented yet");
            }

            @Override
            public double getIntegral(double startTime, double finishTime) {
                throw new UnsupportedOperationException("not implemented yet");
            }

            @Override
            public Type getUnits() {
                return LogisticPopulationSizeModel.this.getUnits();
            }

            @Override
            public void setUnits(Type units) {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Logistic Growth Coalescent";
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

    //
    // protected stuff
    //

    private final Parameter rateParameter;
    private final Parameter shapeParameter;
}
