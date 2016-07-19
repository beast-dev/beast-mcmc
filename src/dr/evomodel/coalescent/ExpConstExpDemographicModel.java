/*
 * ExpConstExpDemographicModel.java
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
import dr.evolution.coalescent.ExpConstExpDemographic;
import dr.evomodelxml.coalescent.ExpConstExpDemographicModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * This class models a two growth-phase demographic with a plateau in the middle
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ExpConstExpDemographicModel.java,v 1.2 2006/08/18 07:44:25 alexei Exp $
 */
public class ExpConstExpDemographicModel extends DemographicModel implements Citable {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ExpConstExpDemographicModel(
            Parameter N0Parameter,
            Parameter N1Parameter,
            Parameter growthRateParameter,
            Parameter timeParameter,
            Parameter relTimeParameter,
            Type units) {

        this(ExpConstExpDemographicModelParser.EXP_CONST_EXP_MODEL, N0Parameter, N1Parameter, growthRateParameter, timeParameter, relTimeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExpConstExpDemographicModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter growthRateParameter, Parameter timeParameter, Parameter relTimeParameter, Type units) {

        super(name);

        expConstExp = new ExpConstExpDemographic(units);

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

        this.relTimeParameter = relTimeParameter;
        addVariable(relTimeParameter);
        relTimeParameter.addBounds(new Parameter.DefaultBounds(1.0, Double.MIN_VALUE, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        expConstExp.setN0(N0Parameter.getParameterValue(0));

        double relTime = relTimeParameter.getParameterValue(0);
        double time2 = timeParameter.getParameterValue(0);

        //System.out.println("relTime=" + relTime);
        //System.out.println("time2=" + (time2));


        double timeInModernGrowthPhase = time2 * relTime;

        double r = -Math.log(N1Parameter.getParameterValue(0)) / timeInModernGrowthPhase;

        //System.out.println("N0=" + N0Parameter.getParameterValue(0));
        //System.out.println("r=" + r);
        //System.out.println("r2=" + growthRateParameter.getParameterValue(0));
        //System.out.println("time1=" + timeInModernGrowthPhase);
        //System.out.println("plateauTime=" + (time2-timeInModernGrowthPhase));

        expConstExp.setGrowthRate(r);
        expConstExp.setGrowthRate2(growthRateParameter.getParameterValue(0));

        expConstExp.setTime1(timeInModernGrowthPhase);
        expConstExp.setPlateauTime(time2 - timeInModernGrowthPhase);

        return expConstExp;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter N1Parameter = null;
    Parameter growthRateParameter = null;
    Parameter timeParameter = null;
    Parameter relTimeParameter = null;
    ExpConstExpDemographic expConstExp = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Exponential-Constant-Exponential Coalescent";
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
}
