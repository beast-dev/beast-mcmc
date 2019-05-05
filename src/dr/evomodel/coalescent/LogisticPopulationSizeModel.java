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

import dr.evomodelxml.coalescent.ExponentialGrowthModelParser;
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
    public LogisticPopulationSizeModel(Parameter N0Parameter, Parameter rateParameter, Parameter shapeParameter, boolean inLogSpace, Type units) {

        this(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, N0Parameter, rateParameter, shapeParameter, inLogSpace, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public LogisticPopulationSizeModel(String name, Parameter N0Parameter, Parameter rateParameter, Parameter shapeParameter, boolean inLogSpace, Type units) {

        super(name, inLogSpace, units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);

        this.rateParameter = rateParameter;
        addVariable(rateParameter);

        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);

        if (isInLogSpace()) {
            N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        } else {
            N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getN0() {
        if (isInLogSpace()) {
            return Math.exp(N0Parameter.getParameterValue(0));
        } else {
            return N0Parameter.getParameterValue(0);
        }
    }

    public double getLogN0() {
        if (isInLogSpace()) {
            return N0Parameter.getParameterValue(0);
        } else {
            return Math.log(N0Parameter.getParameterValue(0));
        }
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

    private final Parameter N0Parameter;
    private final Parameter rateParameter;
    private final Parameter shapeParameter;
}
