package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Package: ALSTreeLikelihood
 * Description:
 *
 *
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 13, 2008
 * Time: 10:13:07 AM
 */
public class ALSTreeLikelihood extends TreeLikelihood
{
    public static final String LIKE_NAME="alsTreeLikelihood";
    protected AbstractObservationProcess observationProcess;

    public ALSTreeLikelihood(AbstractObservationProcess observationProcess, PatternList patternList, TreeModel treeModel, SiteModel siteModel, BranchRateModel branchRateModel, boolean useAmbiguities, boolean storePartials, boolean useScaling) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials, useScaling, false);

        this.observationProcess=observationProcess;
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

    protected double calculateLogLikelihood(){
        // Calculate the partial likelihoods
        super.calculateLogLikelihood();
        // get the frequency model
        double[] freqs=frequencyModel.getFrequencies();
        // let the observationProcess handle the rest
        return observationProcess.nodePatternLikelihood(freqs,branchRateModel,likelihoodCore);
    }


    protected void handleModelChangedEvent(Model model, Object object,  int index){
        if(model==observationProcess || model==treeModel){
                makeDirty();
        }else
            super.handleModelChangedEvent(model,object,index);
    }

    /**
	 * The XML parser
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return LIKE_NAME; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			boolean useAmbiguities = false;
			boolean storePartials = true;
			boolean useScaling = false;
			if (xo.hasAttribute(USE_AMBIGUITIES)) {
				useAmbiguities = xo.getBooleanAttribute(USE_AMBIGUITIES);
			}
			if (xo.hasAttribute(STORE_PARTIALS)) {
				storePartials = xo.getBooleanAttribute(STORE_PARTIALS);
			}
			if (xo.hasAttribute(USE_SCALING)) {
				useScaling = xo.getBooleanAttribute(USE_SCALING);
			}

            AbstractObservationProcess observationProcess= (AbstractObservationProcess)xo.getChild(AbstractObservationProcess.class);
            PatternList patternList = (PatternList)xo.getChild(PatternList.class);
			TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
			SiteModel siteModel = (SiteModel)xo.getChild(SiteModel.class);

            BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);
            Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating ALSTreeLikelihood model.");
            Logger.getLogger("dr.evolution").info("\tIf you publish results using Acquisition-Loss-Mutaion (ALS) Model likelihood, please reference Alekseyenko, Lee and Suchard (in submision).\n---------------------------------\n");

			return new ALSTreeLikelihood(observationProcess, patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials, useScaling);
		}



        //************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of a patternlist on a tree given the site model.";
		}

		public Class getReturnType() { return Likelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
				AttributeRule.newBooleanRule(STORE_PARTIALS, true),
				AttributeRule.newBooleanRule(USE_SCALING, true),
				new ElementRule(PatternList.class),
				new ElementRule(TreeModel.class),
				new ElementRule(SiteModel.class),
				new ElementRule(BranchRateModel.class, true),
                new ElementRule(AbstractObservationProcess.class)
        };
	};

}
