/*
 * MsatFullAncestryImportanceSamplingOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Parameter;
import dr.evolution.tree.NodeRef;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * Produce an importance sample of the ancestry given a msat pattern and a tree.
 */
public class MicrosatelliteFullAncestryImportanceSamplingOperator extends SimpleMCMCOperator {

    public static final String MSAT_FULL_ANCESTRY_IMPORTANCE_SAMPLING_OPERATOR = "MsatFullAncestryImportanceSamplingOperator";
    private Parameter parameter;
    private MicrosatelliteSamplerTreeModel msatSamplerTreeModel;
    private MicrosatelliteModel msatModel;
    private BranchRateModel branchRateModel;


    public MicrosatelliteFullAncestryImportanceSamplingOperator(
            Parameter parameter,
            MicrosatelliteSamplerTreeModel msatSamplerTreeModel,
            MicrosatelliteModel msatModel,
            BranchRateModel branchRateModel,
            double weight){

        super();
        this.parameter = parameter;
        this.msatSamplerTreeModel = msatSamplerTreeModel;
        this.msatModel = msatModel;
        this.branchRateModel = branchRateModel;
        setWeight(weight);
    }

    public double doOperation() {
        TreeModel tree = msatSamplerTreeModel.getTreeModel();

        //get postOrder
        int[] postOrder = new int[tree.getNodeCount()];
        TreeUtils.postOrderTraversalList(tree,postOrder);

        int extNodeCount = tree.getExternalNodeCount();
        double logq=0.0;
        for(int i = 0; i < postOrder.length; i ++){

            //if it's an internal node
            if(postOrder[i] >= extNodeCount){

                //getLikelihoodGiven the children
                NodeRef node = tree.getNode(postOrder[i]);
                NodeRef lc = tree.getChild(node,0);
                NodeRef rc = tree.getChild(node,1);
                int lcState = msatSamplerTreeModel.getNodeValue(lc);
                int rcState = msatSamplerTreeModel.getNodeValue(rc);
                double branchLeftLength = tree.getBranchLength(lc)*branchRateModel.getBranchRate(tree,lc);
                double branchRightLength = tree.getBranchLength(rc)*branchRateModel.getBranchRate(tree,rc);
                double[] probLbranch = msatModel.getColTransitionProbabilities(branchLeftLength, lcState);
                double[] probRbranch = msatModel.getColTransitionProbabilities(branchRightLength, rcState);
                double[] lik = new double[msatModel.getDataType().getStateCount()];
                int currState = (int)parameter.getParameterValue(msatSamplerTreeModel.getParameterIndexFromNodeNumber(postOrder[i]));
                //if node = root node
                if(i == postOrder.length -1){
                    //likelihood of root state also depends on the stationary distribution
                    double[] statDist = msatModel.getStationaryDistribution();
                    for(int j = 0; j < lik.length; j++){
                        lik[j] = probLbranch[j]*probRbranch[j]*statDist[j];
                    }

                }else{

                    for(int j = 0; j < lik.length; j++){
                        lik[j] = probLbranch[j]*probRbranch[j];
                    }

                }

                int sampledState = MathUtils.randomChoicePDF(lik);
                logq = logq + Math.log(lik[currState]) - Math.log(lik[sampledState]);
                parameter.setParameterValue(msatSamplerTreeModel.getParameterIndexFromNodeNumber(postOrder[i]),sampledState);
            }
        }
        
        return logq;
    }

    public String getPerformanceSuggestion(){
        return "None";
    }
    public String getOperatorName(){
        return MSAT_FULL_ANCESTRY_IMPORTANCE_SAMPLING_OPERATOR;
    }
    public int getStepCount(){
        return 1;
    }
}
