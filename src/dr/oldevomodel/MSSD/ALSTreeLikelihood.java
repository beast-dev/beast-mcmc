/*
 * ALSTreeLikelihood.java
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

package dr.oldevomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.LikelihoodScalingProvider;
import dr.oldevomodel.treelikelihood.ScaleFactorsHelper;
import dr.oldevomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;

/**
 * Package: ALSTreeLikelihood
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 13, 2008
 * Time: 10:13:07 AM
 */
public class ALSTreeLikelihood extends TreeLikelihood {

    protected AbstractObservationProcess observationProcess;

    public ALSTreeLikelihood(AbstractObservationProcess observationProcess, PatternList patternList, TreeModel treeModel,
                             SiteModel siteModel, BranchRateModel branchRateModel, boolean useAmbiguities,
                             boolean storePartials, boolean forceRescaling) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials, false, forceRescaling);

        this.observationProcess = observationProcess;
        addModel(observationProcess);

        // TreeLikelihood does not initialize the partials for tips, we'll do it ourselves
        int extNodeCount = treeModel.getExternalNodeCount();
        for (int i = 0; i < extNodeCount; i++) {
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            setPartials(likelihoodCore, patternList, categoryCount, index, i);
        }

/*        //Examine the tree
        double totalTime=0.0;
        double realTime = 0.0;
        for(int i=0; i<treeModel.getNodeCount();++i){
            NodeRef node = treeModel.getNode(i);
            double branchRate = branchRateModel.getBranchRate(treeModel,node);
            double branchTime = treeModel.getBranchLength(node);
            totalTime+=branchRate*branchTime;
            realTime += branchTime;
            System.err.println("Node "+node.toString()+ " time: "+branchTime+ " rate "+branchRate+" together "+branchTime*branchRate);
        }
        System.err.println("TotalTime: "+totalTime);
        System.err.println("RealTime: "+realTime);*/

        scaleFactorsHelper = new ScaleFactorsHelper((LikelihoodScalingProvider) likelihoodCore, this,
                treeModel, stateCount, patternCount, categoryCount);
    }

    protected double calculateLogLikelihood() {
        // Calculate the partial likelihoods
        super.calculateLogLikelihood();
        // get the frequency model
        double[] freqs = frequencyModel.getFrequencies();
        // let the observationProcess handle the rest
        scaleFactorsHelper.resetScaleFactors();
        return observationProcess.nodePatternLikelihood(freqs, likelihoodCore, scaleFactorsHelper);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == observationProcess) {
            likelihoodKnown = false;
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    final private ScaleFactorsHelper scaleFactorsHelper;
}
