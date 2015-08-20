/*
 * TwoEpochDemographicModel.java
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
import dr.evolution.coalescent.TwoEpochDemographic;
import dr.evomodelxml.coalescent.TwoEpochDemographicModelParser;
import dr.inference.model.Parameter;

/**
 * This class models an exponentially growing model that suffers a
 * cataclysmic event and goes into exponential decline
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TwoEpochDemographicModel.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 */
public class TwoEpochDemographicModel extends DemographicModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings.
     */
    public TwoEpochDemographicModel(DemographicModel demo1, DemographicModel demo2, Parameter transitionTimeParameter, Type units) {

        this(TwoEpochDemographicModelParser.TWO_EPOCH_MODEL, demo1, demo2, transitionTimeParameter, units);
    }

    /**
     * Construct demographic model with default settings.
     */
    public TwoEpochDemographicModel(String name, DemographicModel demo1, DemographicModel demo2, Parameter transitionTimeParameter, Type units) {

        super(name);

        this.demo1 = demo1;
        addModel(demo1);
        for (int i = 0; i < demo1.getVariableCount(); i++) {
            addVariable((Parameter) demo1.getVariable(i));
        }

        this.demo2 = demo2;
        addModel(demo2);
        for (int i = 0; i < demo2.getVariableCount(); i++) {
            addVariable((Parameter) demo2.getVariable(i));
        }

        this.transitionTimeParameter = transitionTimeParameter;
        addVariable(transitionTimeParameter);

        setUnits(units);
    }

    // general functions

    public DemographicFunction getDemographicFunction() {

        TwoEpochDemographic twoEpoch = new TwoEpochDemographic(demo1.getDemographicFunction(), demo2.getDemographicFunction(), getUnits());
        twoEpoch.setTransitionTime(transitionTimeParameter.getParameterValue(0));

        return twoEpoch;
    }

    private Parameter transitionTimeParameter = null;
    private DemographicModel demo1 = null;
    private DemographicModel demo2 = null;
}
