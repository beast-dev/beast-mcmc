/*
 * MsatBMA.java
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

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evolution.datatype.Microsatellite;

import java.util.Map;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * This class is used for model averaging over a subset of microsatellite models.
 */
public class MsatBMA extends MicrosatelliteModel{
    private boolean logit;
    public static final int PROP_INDEX = 0;
    public static final int QUAD_INDEX = 1;
    public static final int BIAS_CONST_INDEX = 2;
    public static final int BIAS_LIN_INDEX = 3;
    public static final int GEO_INDEX = 4;
    public static final int PHASE_PROB_INDEX = 5;
    public static final double DEFAULT_VALUE = 0.0;
    public static final int PARAMETER_PRESENT = 1;

    public Parameter[][] paramModelMap;
    public Map<Integer, Integer> modelMap;
    public Parameter modelChoose;
    public Parameter modelIndicator;

    public ArrayList<Parameter> propRates = new ArrayList<Parameter>();
    public ArrayList<Parameter> quadRates = new ArrayList<Parameter>();
    public ArrayList<Parameter> biasConsts = new ArrayList<Parameter>();
    public ArrayList<Parameter> biasLins = new ArrayList<Parameter>();
    public ArrayList<Parameter> geos = new ArrayList<Parameter>();
    public ArrayList<Parameter> phaseProb = new ArrayList<Parameter>();

    public MsatBMA(
            Microsatellite msat,
            boolean logit,
            ArrayList<Parameter> propRates,
            ArrayList<Parameter> quadRates,
            ArrayList<Parameter> biasConsts,
            ArrayList<Parameter> biasLins,
            ArrayList<Parameter> geos,
            ArrayList<Parameter> phaseProb,
            Parameter[][] paramModelMap,
            Parameter modelChoose,
            Parameter modelIndicator,
            Map<Integer, Integer> modelMap){

        super("MsatAveragingModel",msat, null, null);

        for(int i = 0; i < propRates.size(); i++){
            addVariable(propRates.get(i));
        }

        for(int i = 0; i < quadRates.size(); i++){
            addVariable(quadRates.get(i));
        }

        for(int i = 0; i < biasConsts.size(); i++){
            addVariable(biasConsts.get(i));

        }

        for(int i = 0; i < biasLins.size(); i++){
            addVariable(biasLins.get(i));
        }

        for(int i = 0; i < geos.size(); i++){
            addVariable(geos.get(i));
        }

        for(int i = 0; i < phaseProb.size();i++){
            addVariable(phaseProb.get(i));
        }

        addVariable(modelChoose);
        addVariable(modelIndicator);
        this.propRates = propRates;
        this.biasConsts = biasConsts;
        this.biasLins = biasLins;
        this.geos = geos;

        this.logit = logit;
        this.modelChoose = modelChoose;
        this.modelIndicator = modelIndicator;
        this.paramModelMap = paramModelMap;
        this.modelMap = modelMap;

        setupInfinitesimalRates();
        setupMatrix();


    }


    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        modelUpdate = false;
        int paramIndex = -1;
        if(propRates.contains(variable)){
            paramIndex = PROP_INDEX;
        }else if(quadRates.contains(variable)){
            paramIndex = QUAD_INDEX;
        }else if(biasConsts.contains(variable)){
            paramIndex = BIAS_CONST_INDEX;
        }else if(biasLins.contains(variable)){
            paramIndex = BIAS_LIN_INDEX;
        }else if(geos.contains(variable)){
            paramIndex = GEO_INDEX;
        }else if(phaseProb.contains(variable)){
            paramIndex = PHASE_PROB_INDEX;
        }
        //System.out.println("changeModel: "+ variable);
        if(paramIndex > -1){
            //System.out.println(modelMap);
            //System.out.println(getBitVectorValue());
            if(paramModelMap[paramIndex][modelMap.get(getBitVectorValue())] != null){
                updateMatrix = true;
                modelUpdate = true;

            }

        }else if(variable == modelChoose){
            updateMatrix = true;
            indicateModel();
            modelUpdate = true;

        }


    }

    public void indicateModel(){
        modelIndicator.setParameterValueQuietly(0, modelMap.get(getBitVectorValue()));
    }

    public void setupInfinitesimalRates(){
        double rowSum;
        double prop = DEFAULT_VALUE;
        double quad = DEFAULT_VALUE;
        double biasConst = DEFAULT_VALUE;
        double biasLin = DEFAULT_VALUE;
        double geo = DEFAULT_VALUE;
        double phaseProb = DEFAULT_VALUE;

        infinitesimalRateMatrix = new double[stateCount][stateCount];

        //if the rate proportional parameter is present
        if((int)modelChoose.getParameterValue(PROP_INDEX) == PARAMETER_PRESENT){
            //
            prop = getModelParameterValue(PROP_INDEX);
            if((int)modelChoose.getParameterValue(QUAD_INDEX) == PARAMETER_PRESENT){
                quad = getModelParameterValue(QUAD_INDEX);
            }
        }
        for(int i = 0; i < stateCount;i++){
            rowSum = 0.0;
            if(i - 1 > -1){
                infinitesimalRateMatrix[i][i - 1] = 1+prop*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i - 1];
            }

            if(i + 1 < stateCount){
                infinitesimalRateMatrix[i][i + 1] = 1+prop*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i + 1];
            }
            infinitesimalRateMatrix[i][i] = -rowSum;

        }

        if((int)modelChoose.getParameterValue(BIAS_CONST_INDEX) == PARAMETER_PRESENT){

            biasConst = getModelParameterValue(BIAS_CONST_INDEX);
            //System.out.print("biasConst: "+biasConst);
            if((int)modelChoose.getParameterValue(BIAS_LIN_INDEX) == PARAMETER_PRESENT){
                biasLin = getModelParameterValue(BIAS_LIN_INDEX);
                //System.out.print("biasLin: "+biasLin);
            }

            double[][] subRates =  infinitesimalRateMatrix;
            infinitesimalRateMatrix = new double[stateCount][stateCount];

            LinearBiasModel.setupInfinitesimalRates(
                    infinitesimalRateMatrix,
                    subRates,
                    biasConst,
                    biasLin,
                    stateCount,
                    logit
            );
        }

        if((int)modelChoose.getParameterValue(GEO_INDEX) == PARAMETER_PRESENT){
            geo = getModelParameterValue(GEO_INDEX);
            if((int)modelChoose.getParameterValue(PHASE_PROB_INDEX) == PARAMETER_PRESENT){
                phaseProb = getModelParameterValue(PHASE_PROB_INDEX);
            }
            double[][] subRates =  infinitesimalRateMatrix;
            infinitesimalRateMatrix = new double[stateCount][stateCount];

            TwoPhaseModel.setupInfinitesimalRates(
                stateCount,
                geo,
                phaseProb,
                infinitesimalRateMatrix,
                subRates
            );
        }
    }

    private double getModelParameterValue(int paramIndex){
        int modelCode = modelMap.get(getBitVectorValue());
        return (paramModelMap[paramIndex][modelCode]).getParameterValue(0);

    }

    private int getBitVectorValue(){
        String bitVec = "";
        for(int i = 0; i < modelChoose.getDimension(); i++){
            bitVec = bitVec + (int)modelChoose.getParameterValue(i);
        }
        

        return Integer.parseInt(bitVec,2);
    }

    public void computeStationaryDistribution(){
        if((int)modelChoose.getParameterValue(GEO_INDEX) == PARAMETER_PRESENT){
            computeTwoPhaseStationaryDistribution();
        }else{
            computeOnePhaseStationaryDistribution();
        }
        super.computeStationaryDistribution();
    }

}
