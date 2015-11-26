/*
 * ScaledPiecewiseModel.java
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
import dr.evolution.coalescent.PiecewiseLinearPopulation;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ScaledPiecewiseModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ScaledPiecewiseModel.java,v 1.3 2005/04/11 11:24:39 alexei Exp $
 */
public class ScaledPiecewiseModel extends DemographicModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ScaledPiecewiseModel(Parameter N0Parameter, TreeModel treeModel, boolean linear, Type units) {

        this(ScaledPiecewiseModelParser.PIECEWISE_POPULATION, N0Parameter, treeModel, linear, units);

    }

    /**
     * Construct demographic model with default settings
     */
    public ScaledPiecewiseModel(String name, Parameter N0Parameter, TreeModel treeModel, boolean linear, Type units) {

        super(name);

        if (N0Parameter.getDimension() < 2) {
            throw new IllegalArgumentException("Must have at least 2 epochs");
        }

        this.N0Parameter = N0Parameter;
        this.treeModel = treeModel;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, N0Parameter.getDimension()));

        setUnits(units);

        int epochs = N0Parameter.getDimension();

        if (linear) {
            piecewiseFunction = new PiecewiseLinearPopulation(new double[epochs - 1], new double[epochs], units);
        } else {
            piecewiseFunction = new PiecewiseConstantPopulation(new double[epochs - 1], new double[epochs], units);
        }
    }

    /**
     *
     */
    public DemographicFunction getDemographicFunction() {

        double height = treeModel.getNodeHeight(treeModel.getRoot());
        int epochs = N0Parameter.getDimension();

        for (int i = 0; i < N0Parameter.getDimension(); i++) {
            piecewiseFunction.setArgument(i, N0Parameter.getParameterValue(i));
            if (i < epochs - 1) piecewiseFunction.setEpochDuration(i, height / (double) epochs);
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

    //
    // protected stuff
    //

    Parameter N0Parameter;
    TreeModel treeModel;
    PiecewiseConstantPopulation piecewiseFunction = null;
}
