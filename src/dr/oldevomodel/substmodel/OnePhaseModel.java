/*
 * OnePhaseModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;


/**
 * @author Chieh-Hsi Wu
 *
 * An abstract for One Phase Microsatellite Models
 */
public abstract class OnePhaseModel extends MicrosatelliteModel{
    protected ArrayList<Variable<Double>> nestedParams = null;

    /**
     * Constructor
     *
     * @param name              Model name
     * @param microsatellite    Microsatellite data type
     * @param freqModel         Equilibrium frequencies
     * @param parameter         Infinitesimal rates
     */
    public OnePhaseModel(String name, Microsatellite microsatellite, FrequencyModel freqModel, Parameter parameter){
        super(name, microsatellite, freqModel, parameter);
        nestedParams=new ArrayList<Variable<Double>>();
    }

    /*
     * adding the parameters only if its not a submodel.
     */
    protected void addParam(Variable<Double> param){
        if(isNested){
            nestedParams.add(param);
        }else{
            super.addVariable(param);
        }
    }

    /*
     * get the parameters in this submodel
     */
    public Variable<Double> getNestedParameter(int i){
        return nestedParams.get(i);
    }

    /*
     * get number of nested parameters in the submodel
     */
    public int getNestedParameterCount(){
        return nestedParams.size();
    }

    /*
     * The One Phase Models are special cases of the birth-death chain,
     * and therefore we can use this to calculate the stationay distribution
     * given a infinitesimal rate matrix.
     */
    public void computeStationaryDistribution(){
        if(useStationaryFreqs){
            computeOnePhaseStationaryDistribution();
        }
        super.computeStationaryDistribution();

    }


}
