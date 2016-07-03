/*
 * LinearBiasModel.java
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
 * Implements the LinearBiasModel of microsatellites.
 */
public class LinearBiasModel extends OnePhaseModel{

    /*
     *Parameters for setting up the infinitesimal rate matrix.
     */
    private Parameter biasConst;
    private Parameter biasLin;
    private ArrayList<Variable<Double>> submodelParameters = null;
    private boolean estimateSubmodelParams = false;
    private boolean updateSubmodelRates = false;
    private boolean inLogitSpace = false;

    public static final double delta = 1e-15;

    public static final String LINEAR_BIAS_MODEL = "LINEARBIASModel";


    /**
     * Constructor
     *
     * @param microsatellite            microsatellite data type
     * @param freqModel                 user specified initial equilibrium frequencies
     * @param submodel                  submodel of this linear bias model
     * @param biasConst                 bias constant parameter
     * @param biasLinear                bias linear parameter
     * @param inLogitSpace              indicates whether the bias parameters are in the logit space
     * @param estimateSubmodelParams    inidicate whether the parameters of submodel will be estimated
     * @param isSubmodel                inidicate whether this model is a submodel of another microsatellite model
     */
    public LinearBiasModel(
            Microsatellite microsatellite,
            FrequencyModel freqModel,
            OnePhaseModel submodel,
            Parameter biasConst,
            Parameter biasLinear,
            boolean inLogitSpace,
            boolean estimateSubmodelParams,
            boolean isSubmodel){

        super(LINEAR_BIAS_MODEL, microsatellite, freqModel, null);

        isNested = isSubmodel;
        this.subModel = submodel;
        this.estimateSubmodelParams = estimateSubmodelParams;
        if(this.estimateSubmodelParams){
            submodelParameters = new ArrayList<Variable<Double>>();
            for(int i = 0; i < subModel.getNestedParameterCount(); i++){

                if(isNested){
                    addVariable(subModel.getNestedParameter(i));
                }

                addParam(subModel.getNestedParameter(i));
                submodelParameters.add(subModel.getNestedParameter(i));
            }
            updateSubmodelRates = true;
        }



        //The default setting of the parameters gives infinitesimal rates with no directional bias.
        if(biasConst != null){
            this.biasConst = biasConst;
        }else{
            this.biasConst = new Parameter.Default(0.5);
        }


        if(biasLinear != null){
            biasLin = biasLinear;
        }else{
            biasLin = new Parameter.Default(0.0);
        }
        addParam(this.biasConst);
        addParam(this.biasLin);

        this.inLogitSpace = inLogitSpace;

        //printDetails();

        setupInfinitesimalRates();

        if(freqModel == null){
            useStationaryFreqs = true;
            computeStationaryDistribution();
        }else{
            this.freqModel = freqModel;

        }
        addModel(this.freqModel);

    }
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(submodelParameters !=null && submodelParameters.indexOf((Parameter)variable) != -1){
            updateSubmodelRates = true;
        }
        updateMatrix = true;

    }

     public void setupInfinitesimalRates(){
        if(updateSubmodelRates){
            subModel.setupInfinitesimalRates();
            updateSubmodelRates = false;
        }

        double biasConst = this.biasConst.getParameterValue(0);
        double biasLin = this.biasLin.getParameterValue(0);
        setupInfinitesimalRates(
                infinitesimalRateMatrix,
                subModel.getInfinitesimalRates(),
                biasConst,
                biasLin,
                stateCount,
                inLogitSpace
        );

    }


    public static void setupInfinitesimalRates(
            double[][] rates,
            double[][] subModelRateMatrix,
            double biasConst,
            double biasLin,
            int stateCount,
            boolean inLogitSpace){
        double rowSum;
        double expansionProb = 0.5;
        for(int i = 0; i < stateCount;i++){
            rowSum = 0.0;
            expansionProb = computeExpansionProb(biasConst,biasLin,i, inLogitSpace);
            if(expansionProb < delta){
                System.out.println("changing expan prob from " + expansionProb+ " to " + delta
                +"\nbiasConst: "+biasConst+", biasLin: "+biasLin);
                expansionProb = delta;

            }else if (expansionProb > (1.0-delta)){

                System.out.println("changing expan prob from " + expansionProb+ " to " + (1.0-delta)
                        +"\nbiasConst: "+biasConst+", biasLin: "+biasLin);
                expansionProb = 1.0-delta;
            }
            if(i - 1 > -1){
                rates[i][i - 1] = subModelRateMatrix[i][i-1]*(1.0 - expansionProb);
                rowSum = rowSum+rates[i][i - 1];
            }
            if(i + 1 < stateCount){
                rates[i][i + 1] = subModelRateMatrix[i][i+1]*expansionProb;
                rowSum = rowSum + rates[i][i + 1];
            }

            rates[i][i] = rowSum*-1;

        }
    }

    public static double computeExpansionProb(double biasConst, double biasLin, int length, boolean inLogitSpace){
        double expanProb = 0.5;
        if(inLogitSpace){
            double numerator = Math.exp(biasConst+biasLin*length);
            expanProb = numerator/(1+numerator);
        }else{
            expanProb = biasConst+biasLin*length;
        }
        return  expanProb;

    }


    public Parameter getBiasConstant(){
        return biasConst;
    }

    public Parameter getBiasLinearPercent(){
        return biasLin;
    }

    public boolean isEstimatingSubmodelParams(){
        return estimateSubmodelParams;
    }

    public boolean isInLogitSpace(){
        return inLogitSpace;
    }

    public void printDetails(){
        System.out.println("Details of the Linear Bias Model and its paramters:");
        System.out.println("a submodel:                     "+isNested);
        System.out.println("in logit space:                 "+inLogitSpace);
        System.out.println("has submodel:                   "+hasSubmodel());
        if(hasSubmodel()){
            System.out.println("submodel class:                 "+subModel.getClass());
        }
        System.out.println("esitmating submodel parameters: "+estimateSubmodelParams);
        System.out.println("bias constant:                  "+biasConst.getParameterValue(0));
        System.out.println("bias linear coefficient:        "+biasLin.getParameterValue(0));
    }

}
