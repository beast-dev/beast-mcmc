/*
 * MsatSingleAncestralStateGibbsOperator.java
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

import dr.inference.model.Parameter;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evolution.tree.NodeRef;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 *
 */
public class MicrosatelliteSingleAncestralStateGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    public static final String MSAT_SINGLE_ANCESTAL_STATE_GIBBS_OPERATOR = "MsatSingleAncestralStateGibbsOperator";
    private Parameter parameter;
    private MicrosatelliteSamplerTreeModel msatSamplerTreeModel;
    private MicrosatelliteModel msatModel;
    private BranchRateModel branchRateModel;

    public MicrosatelliteSingleAncestralStateGibbsOperator(
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

        
        int index = MathUtils.nextInt(parameter.getDimension());

        //double logq=0.0;


        //getLikelihoodGiven the children
        NodeRef node = tree.getNode(msatSamplerTreeModel.getNodeNumberFromParameterIndex(index));
        NodeRef lc = tree.getChild(node,0);
        NodeRef rc = tree.getChild(node,1);


        int lcState = msatSamplerTreeModel.getNodeValue(lc);
        int rcState = msatSamplerTreeModel.getNodeValue(rc);


        double branchLeftLength = tree.getBranchLength(lc)*branchRateModel.getBranchRate(tree,lc);
        double branchRightLength = tree.getBranchLength(rc)*branchRateModel.getBranchRate(tree,rc);


        double[] probLbranch = msatModel.getColTransitionProbabilities(branchLeftLength, lcState);
        double[] probRbranch = msatModel.getColTransitionProbabilities(branchRightLength, rcState);

        double[] lik = new double[msatModel.getDataType().getStateCount()];
        //int currState = (int)parameter.getParameterValue(index);
        //if node = root node

        if(tree.isRoot(node)){
            //likelihood of root state also depends on the stationary distribution
            double[] statDist = msatModel.getStationaryDistribution();
            for(int j = 0; j < lik.length; j++){
                lik[j] = probLbranch[j]*probRbranch[j]*statDist[j];
            }
        }else{
            NodeRef parent = tree.getParent(node);
            int pState = msatSamplerTreeModel.getNodeValue(parent);
            double branchParentLength = tree.getBranchLength(node)*branchRateModel.getBranchRate(tree,node);
            double[] probPbranch = msatModel.getRowTransitionProbabilities(branchParentLength,pState);
            for(int j = 0; j < lik.length; j++){
                lik[j] = probLbranch[j]*probRbranch[j]*probPbranch[j];
            }
        }

        int sampledState = MathUtils.randomChoicePDF(lik);
        //logq = logq + Math.log(lik[currState]) - Math.log(lik[sampledState]);
        parameter.setParameterValue(index,sampledState);
        return 0.0;
    }

    public String getPerformanceSuggestion(){
        return "None";
    }
    public String getOperatorName(){
        return MSAT_SINGLE_ANCESTAL_STATE_GIBBS_OPERATOR;
    }
    public int getStepCount(){
        return 1;
    }

}