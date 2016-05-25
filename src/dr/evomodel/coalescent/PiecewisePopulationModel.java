/*
 * PiecewisePopulationModel.java
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
import dr.evolution.coalescent.PiecewiseConstantPopulation;
import dr.evolution.coalescent.PiecewiseExponentialPopulation;
import dr.evolution.coalescent.PiecewiseLinearPopulation;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: PiecewisePopulationModel.java,v 1.13 2005/05/24 20:25:57 rambaut Exp $
 */
public class PiecewisePopulationModel extends DemographicModel {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public PiecewisePopulationModel(String name, Parameter N0Parameter, double[] epochLengths, boolean isLinear, Type units) {

        super(name);

        this.epochCount = epochLengths.length + 1;

        if (N0Parameter.getDimension() != epochCount) {
            throw new IllegalArgumentException(
                    "epochSize parameter must have the same dimensions as the number of epochs: (" + epochCount +
                            ") but instead has " + N0Parameter.getDimension() + "!"
            );
        }

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, N0Parameter.getDimension()));

        //addVariable(epochLengths);
        //epochLengths.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, epochLengths.getDimension()));

        setUnits(units);

        if (isLinear) {
            piecewiseFunction = new PiecewiseLinearPopulation(epochLengths, new double[N0Parameter.getDimension()], units);
        } else {
            piecewiseFunction = new PiecewiseConstantPopulation(epochLengths, new double[N0Parameter.getDimension()], units);
        }
    }

    /**
     * Construct demographic model with default settings
     */
    public PiecewisePopulationModel(String name, Parameter N0Parameter, Parameter growthRatesParameter,
                                    double[] epochLengths, Type units) {

        super(name);

        this.epochCount = epochLengths.length + 1;

        this.N0Parameter = N0Parameter;
        this.growthRatesParameter = growthRatesParameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, N0Parameter.getDimension()));

        addVariable(growthRatesParameter);
        growthRatesParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, growthRatesParameter.getDimension()));

        setUnits(units);

        int popSizeCount = N0Parameter.getDimension();
        int growthRateCount = growthRatesParameter.getDimension();

        if (popSizeCount == epochCount && growthRateCount == 1) {
            piecewiseFunction = new PiecewiseExponentialPopulation(epochLengths,
                    new double[N0Parameter.getDimension()],
                    growthRatesParameter.getParameterValue(0),
                    units);
        } else if (popSizeCount == 1 && growthRateCount == epochCount) {
            piecewiseFunction = new PiecewiseExponentialPopulation(epochLengths,
                    N0Parameter.getParameterValue(0),
                    new double[growthRatesParameter.getDimension()], units);
        } else {
            if (growthRatesParameter.getDimension() != epochCount) {
                throw new IllegalArgumentException(
                        "growthRate parameter must have the same dimension as the number of epochs: (" + epochCount +
                                ") but instead has " + N0Parameter.getDimension() + "!"
                );
            }
        }


        addStatistic(new GrowthRateStatistic());
    }

    public DemographicFunction getDemographicFunction() {
        if (growthRatesParameter != null) {
            // exponential growth
            for (int i = 0; i < N0Parameter.getDimension(); i++) {
                piecewiseFunction.setArgument(i, N0Parameter.getParameterValue(i));
            }
            for (int i = 0; i < growthRatesParameter.getDimension(); i++) {
                piecewiseFunction.setArgument(i + N0Parameter.getDimension(), growthRatesParameter.getParameterValue(i));
            }
        } else {
            // constant or linear growth
            for (int i = 0; i < N0Parameter.getDimension(); i++) {
                piecewiseFunction.setArgument(i, N0Parameter.getParameterValue(i));
            }
        }
        return piecewiseFunction;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {

        if (parameter == N0Parameter) {
            //System.out.println("popSize parameter changed..");
        }

        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    public class GrowthRateStatistic extends Statistic.Abstract {

        public GrowthRateStatistic() {
            super("growthRate");
        }

        public int getDimension() {
            return ((PiecewiseExponentialPopulation) piecewiseFunction).getEpochCount();
        }

        public double getStatisticValue(int i) {
            return ((PiecewiseExponentialPopulation) piecewiseFunction).getEpochGrowthRate(i);
        }

    }

    //
    // private stuff
    //

    private Parameter N0Parameter;
    private Parameter growthRatesParameter;
    private DemographicFunction piecewiseFunction = null;

    private final int epochCount;
}
