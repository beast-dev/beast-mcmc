/*
 * CataclysmicDemographicModel.java
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

import dr.evolution.coalescent.CataclysmicDemographic;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodelxml.coalescent.CataclysmicDemographicModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * This class models an exponentially growing model that suffers a
 * cataclysmic event and goes into exponential decline
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: CataclysmicDemographicModel.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 */
public class CataclysmicDemographicModel extends DemographicModel implements Citable {

    /**
     * Construct demographic model with default settings
     */
    public CataclysmicDemographicModel(Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter,
                Type units, boolean useSpike) {

        this(CataclysmicDemographicModelParser.CATACLYSM_MODEL, N0Parameter, N1Parameter, growthRateParameter, timeParameter, units, useSpike);
    }

    /**
     * Construct demographic model with default settings
     */
    public CataclysmicDemographicModel(String name, Parameter N0Parameter, Parameter secondParam, Parameter growthRateParameter,
                Parameter timeParameter, Type units, boolean useSpike) {

        super(name);

        cataclysm = new CataclysmicDemographic(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        if (useSpike) {
            this.N1Parameter = secondParam;
            N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            addVariable(N1Parameter);
        } else {
            this.declineRateParameter = secondParam;
            addVariable(declineRateParameter);
        }

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
        double declineRate = (declineRateParameter == null) ? Math.log(N1Parameter.getParameterValue(0)) / t :
                declineRateParameter.getParameterValue(0);
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
    Parameter declineRateParameter = null;
    CataclysmicDemographic cataclysm = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Boom-Bust Coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("B", "Shapiro"),
                    new Author("", "et al.")
            },
            "Rise and fall of the Beringian steppe bison",
            2004,
            "Science",
            306, 1561, 1565
    );
}
