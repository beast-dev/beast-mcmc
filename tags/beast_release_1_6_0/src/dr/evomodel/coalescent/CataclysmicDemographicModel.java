/*
 * CataclysmicDemographicModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.CataclysmicDemographic;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodelxml.coalescent.CataclysmicDemographicModelParser;
import dr.inference.model.Parameter;

/**
 * This class models an exponentially growing model that suffers a
 * cataclysmic event and goes into exponential decline
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: CataclysmicDemographicModel.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 */
public class CataclysmicDemographicModel extends DemographicModel {

    /**
     * Construct demographic model with default settings
     */
    public CataclysmicDemographicModel(Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter, Type units) {

        this(CataclysmicDemographicModelParser.CATACLYSM_MODEL, N0Parameter, N1Parameter, growthRateParameter, timeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public CataclysmicDemographicModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter, Type units) {

        super(name);

        cataclysm = new CataclysmicDemographic(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.N1Parameter = N1Parameter;
        addVariable(N1Parameter);
        N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));

        this.timeParameter = timeParameter;
        addVariable(timeParameter);
        timeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));


        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        cataclysm.setN0(N0Parameter.getParameterValue(0));
        cataclysm.setGrowthRate(growthRateParameter.getParameterValue(0));
        cataclysm.setCataclysmTime(timeParameter.getParameterValue(0));

        // Doesn't this...
        /*
          double N0 = N0Parameter.getParameterValue(0);
          double N1 = N1Parameter.getParameterValue(0) * N0;
          double t = timeParameter.getParameterValue(0);
          double declineRate = Math.log(N1/N0)/t;
          */ // ..collapse to...

        double t = timeParameter.getParameterValue(0);
        double declineRate = Math.log(N1Parameter.getParameterValue(0)) / t;
        cataclysm.setDeclineRate(declineRate);


        return cataclysm;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter N1Parameter = null;
    Parameter growthRateParameter = null;
    Parameter timeParameter = null;
    CataclysmicDemographic cataclysm = null;
}
