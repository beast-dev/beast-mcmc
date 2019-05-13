/*
 * DemographicModel.java
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

import dr.evolution.util.Units;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * This interface provides a model that returns a PopulationSizeFunction.
 *
 * @author Andrew Rambaut
 */
public abstract class PopulationSizeModel extends AbstractModel implements Units {

    public PopulationSizeModel(String name, Parameter logN0Parameter, Type units) {
        super(name);

        if (logN0Parameter != null) {
            this.logN0Parameter = logN0Parameter;
            logN0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        } else {
            // used to store the parameter value when used in a multi-epoch chain
            this.logN0Parameter = new Parameter.Default(1);
        }
        addVariable(this.logN0Parameter);

        setUnits(units);
    }

    // general functions

    public abstract PopulationSizeFunction getPopulationSizeFunction();

    /**
     * Package private because this is used by the PiecewisePopulationSizeModel to chain epochs
     * together. Sets the parameter quietly to avoid sending further update messages.
     *
     * @param logN0 the log N0
     */
    final void setLogN0(double logN0) {
        logN0Parameter.setParameterValueQuietly(0, logN0);

    }

    public final double getN0() {
        return Math.exp(logN0Parameter.getParameterValue(0));
    }

    public final double getLogN0() {
        return logN0Parameter.getParameterValue(0);
    }

    private final Parameter logN0Parameter;

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Units in which population size is measured.
     */
    private Type units;

    /**
     * sets units of measurement.
     *
     * @param u units
     */
    public void setUnits(Type u) {
        units = u;
    }

    /**
     * returns units of measurement.
     */
    public Type getUnits() {
        return units;
    }

}