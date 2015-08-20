/*
 * EmergingEpidemicModel.java
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
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.EmergingEpidemicModelParser;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

/**
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ExponentialGrowthModel.java,v 1.14 2005/05/24 20:25:57 rambaut Exp $
 */
public class EmergingEpidemicModel extends DemographicModel {

    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public EmergingEpidemicModel(Parameter growthRateParameter,
                                 Parameter generationTimeParameter,
                                 Parameter generationShapeParameter,
                                 Parameter offspringDispersionParameter,
                                 TreeModel treeModel,
                                 Type units) {

        this(EmergingEpidemicModelParser.EMERGING_EPIDEMIC_MODEL, growthRateParameter, generationTimeParameter, generationShapeParameter, offspringDispersionParameter, treeModel, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public EmergingEpidemicModel(String name,
                                 Parameter growthRateParameter,
                                 Parameter generationTimeParameter,
                                 Parameter generationShapeParameter,
                                 Parameter offspringDispersionParameter,
                                 TreeModel treeModel,
                                 Type units) {

        super(name);

        exponentialGrowth = new ExponentialGrowth(units);

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.generationTimeParameter = generationTimeParameter;
        addVariable(generationTimeParameter);
        generationTimeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.generationShapeParameter = generationShapeParameter;
        addVariable(generationShapeParameter);
        generationShapeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.offspringDispersionParameter = offspringDispersionParameter;
        addVariable(offspringDispersionParameter);
        offspringDispersionParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.treeModel = treeModel;
        addModel(treeModel);

        addStatistic(new N0Statistic("N0"));
        addStatistic(new RStatistic("R"));

        setUnits(units);
    }

    // general functions

    public DemographicFunction getDemographicFunction() {
        exponentialGrowth.setN0(getN0());
        exponentialGrowth.setGrowthRate(growthRateParameter.getParameterValue(0));

        return exponentialGrowth;
    }

    public double getR() {
        double r = growthRateParameter.getParameterValue(0);
        double Tg = generationTimeParameter.getParameterValue(0);
        double alpha = generationShapeParameter.getParameterValue(0);

        double R = Math.pow(1.0 + ((r * Tg) / alpha), alpha);

        return R;
    }

    public double getN0() {
        double R = getR();

        double t0 = treeModel.getNodeHeight(treeModel.getRoot());

        double r = growthRateParameter.getParameterValue(0);
        double Tg = generationTimeParameter.getParameterValue(0);
        double k = offspringDispersionParameter.getParameterValue(0);

        double N0 = (k * Tg * Math.exp(r * t0)) / (R * (k + R));

        return N0;
    }

    public class N0Statistic extends Statistic.Abstract {

        public N0Statistic(String name) {
            super(name);
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(final int i) {
            return getN0();
        }
    }

    public class RStatistic extends Statistic.Abstract {

        public RStatistic(String name) {
            super(name);
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(final int i) {
            return getR();
        }
    }

    //
    // protected stuff
    //

    private final Parameter growthRateParameter;
    private final Parameter generationTimeParameter;
    private final Parameter generationShapeParameter;
    private final Parameter offspringDispersionParameter;
    private final TreeModel treeModel;

    private final ExponentialGrowth exponentialGrowth;
}
