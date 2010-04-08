package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
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
                             SiteModel siteModel, BranchRateModel branchRateModel, boolean useAmbiguities, boolean storePartials) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials, false, false);

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

    }

    protected double calculateLogLikelihood() {
        // Calculate the partial likelihoods
        super.calculateLogLikelihood();
        // get the frequency model
        double[] freqs = frequencyModel.getFrequencies();
        // let the observationProcess handle the rest
        return observationProcess.nodePatternLikelihood(freqs, likelihoodCore);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == observationProcess) {
            likelihoodKnown = false;            
        } else
            super.handleModelChangedEvent(model, object, index);
    }
}
