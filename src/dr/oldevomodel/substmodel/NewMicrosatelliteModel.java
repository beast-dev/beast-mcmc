/*
 * NewMicrosatelliteModel.java
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
import dr.evolution.datatype.Microsatellite;
import dr.math.ModifiedBesselFirstKind;

/**
 * @author Chieh-Hsi Wu
 * Implementation of models by Watkins (2007)
 */
public class NewMicrosatelliteModel extends MicrosatelliteModel {

    Parameter biasConst;
    private boolean normalize;
    public NewMicrosatelliteModel(Microsatellite msat, FrequencyModel rootFreqModel){
        this(msat, rootFreqModel, false);
        biasConst = new Parameter.Default(0.5);
    }
    public NewMicrosatelliteModel(Microsatellite msat, FrequencyModel rootFreqModel, boolean normalize){
        super("NewMicrosatelliteModel", msat, rootFreqModel,null);
        this.normalize = normalize;
        biasConst = new Parameter.Default(0.5);
        double[] stationaryDist = new double[stateCount];
        for(int i = 0; i < stationaryDist.length;i++){
            stationaryDist[i] = 1.0/stateCount;
        }
        freqModel = new FrequencyModel(dataType, stationaryDist);
        computeStationaryDistribution();
    }
    protected void storeState(){};
    protected void restoreState(){};
    public void getTransitionProbabilities(double distance, double[] matrix){
        int k = 0;
        double[] rowSums = new double[stateCount];
        double bCVal = biasConst.getParameterValue(0);
        for(int i = 0; i < stateCount; i ++){
            for(int j = 0; j < stateCount; j++){
                int n = i - j;
                //matrix[k] = Math.exp(-distance)* Math.pow(bCVal/(1-bCVal),n/2.0)*ModifiedBesselFirstKind.bessi(2*Math.sqrt(bCVal*(1-bCVal))*distance,Math.abs(n));
                matrix[k] = Math.exp(-distance)*ModifiedBesselFirstKind.bessi(distance,Math.abs(n));
                rowSums[i] += matrix[k];
                k++;
            }
            //System.out.println(rowSums[i]);

        }
        if(normalize){
            k = 0;
            for(int i = 0; i < stateCount; i ++){
                for(int j = 0; j < stateCount; j++){
                    matrix[k] =  matrix[k]/rowSums[i];
                    k++;
                }
            }
        }
    }

    public double[] getRowTransitionProbabilities(double distance, int parentState){

        double[] probabilities = new double[stateCount];
        for(int i = 0; i < probabilities.length;i++){
            int n = parentState - i;
            probabilities[i] = Math.exp(-distance)*ModifiedBesselFirstKind.bessi(distance,Math.abs(n));
        }

        return probabilities;
    }
   
    public double[] getColTransitionProbabilities(double distance, int childState){
        double[] probabilities = new double[stateCount];
        for(int i = 0; i < probabilities.length;i++){
            int n = i - childState;
            probabilities[i] = Math.exp(-distance)*ModifiedBesselFirstKind.bessi(distance,Math.abs(n));
        }
        return probabilities;
    }
    public double getLogOneTransitionProbabilityEntry(double distance, int parentState, int childState){
        return Math.log(getOneTransitionProbabilityEntry(distance, parentState, childState));
    }

    public double getOneTransitionProbabilityEntry(double distance, int parentState, int childState){
        int n = parentState - childState;
        double probability = Math.exp(-distance)*ModifiedBesselFirstKind.bessi(distance,Math.abs(n));
        return probability;
    }
    protected void ratesChanged() {};
    protected void setupRelativeRates(){};
    public void setupInfinitesimalRates(){};
    protected void frequenciesChanged() {};

    public static void main(String[] args){
        Microsatellite msat = new Microsatellite(1,5);
        NewMicrosatelliteModel nmsatModel = new NewMicrosatelliteModel(msat, null);
        double[] probs = new double[msat.getStateCount()*msat.getStateCount()];
        nmsatModel.getTransitionProbabilities(1.0,probs);
        int k =0;
        for(int i = 0; i < msat.getStateCount(); i++){
            for(int j = 0; j < msat.getStateCount(); j++){
                System.out.print(probs[k++]+" ");
            }
            System.out.println();
        }
        double[] statDist = nmsatModel.getStationaryDistribution();
        for(int i = 0; i < statDist.length; i++){
            System.out.print(statDist[i]+" ");
        }
    }

}
