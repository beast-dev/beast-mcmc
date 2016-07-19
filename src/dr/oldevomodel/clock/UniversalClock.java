/*
 * UniversalClock.java
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

package dr.oldevomodel.clock;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * A class that calculates the rate of evolution based on of the mass and temperature
 *
 * @author Alexei Drummond
 * @version $Id: UniversalClock.java,v 1.5 2005/02/21 15:16:09 rambaut Exp $
 */

public class UniversalClock extends AbstractBranchRateModel {

    public static final String UNIVERSAL_CLOCK = "universalClock";

    /**
     * creates a universal clock model
     *
     * @param rateParameter        the rate vector of *all* nodes
     * @param massParameter        the mass vector of *all* nodes
     * @param temperatureParameter the temperature vector of *all* nodes
     * @param scaleParameter       a single-dimensional scale parameter
     */
    public UniversalClock(
            Parameter rateParameter,
            Parameter massParameter,
            Parameter temperatureParameter,
            Parameter scaleParameter
    ) {

        super(UNIVERSAL_CLOCK);

        this.rateParameter = rateParameter;
        this.massParameter = massParameter;
        this.temperatureParameter = temperatureParameter;
        this.scaleParameter = scaleParameter;

        // don't add rate parameter, cause that is what you are changing!
        // you don't care if it changes
        addVariable(massParameter);
        addVariable(temperatureParameter);
        addVariable(scaleParameter);

    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no submodels so nothing to do
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        if ((variable == massParameter) || (variable == temperatureParameter)) {
            if (index == -1) {
                calculateAllRates();
            } else calculateRate(index);
        } else if (variable == scaleParameter) {
            calculateAllRates();
        } else {
            throw new RuntimeException("unknown parameter changed in " + UNIVERSAL_CLOCK);
        }
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // Private methods
    // **************************************************************

    private void calculateAllRates() {
        int numNodes = massParameter.getDimension();
        for (int i = 0; i < numNodes; i++) {
            calculateRate(i);
        }
    }

    private void calculateRate(int index) {

        double mass = massParameter.getParameterValue(index);
        double temperature = temperatureParameter.getParameterValue(index);

        double scale = scaleParameter.getParameterValue(0);

        // replace this with the real equation!!
        double substitutionRate = scale * Math.pow(mass, -0.25) * Math.exp(temperature);

        rateParameter.setParameterValue(index, substitutionRate);
    }

    public double getBranchRate(Tree tree, NodeRef node) {
        throw new RuntimeException("Look at code before running this class!");
    }

    Parameter rateParameter = null;
    Parameter massParameter = null;
    Parameter temperatureParameter = null;
    Parameter scaleParameter = null;

}

