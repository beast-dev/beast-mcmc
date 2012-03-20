package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

public class BeagleSequenceSimulator {

	BeagleSequenceSimulator(TreeModel treeModel, //
			BranchSubstitutionModel branchSubstitutionModel, //
			SiteRateModel siteRateModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int sequenceLength //
	) {

		  int tipCount = treeModel.getExternalNodeCount();
          int nodeCount = treeModel.getNodeCount();
          int eigenCount = branchSubstitutionModel.getEigenCount();
          int internalNodeCount = treeModel.getInternalNodeCount();
          int scaleBufferCount = internalNodeCount + 1;
		
        // TODO: useAmbiguities ?
        int compactPartialsCount = tipCount;
        
        // TODO: do we need this when simulating?
        int patternCount = sequenceLength;
        
        // TODO: don't really need a siteModel to get the stateCount
        // too bad the siteModel only functions as a container with respect to the substitutionModel
        // parse from BranchSubstitutionModel ?
        int stateCount = freqModel.getDataType().getStateCount();
        
        // TODO: don't like this categoryCount though, does require a siteModel
        int categoryCount = siteRateModel.getCategoryCount();
        
        // one partials buffer for each tip and two for each internal node (for store restore)
        BufferIndexHelper partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

        // two eigen buffers for each decomposition for store and restore.
        BufferIndexHelper eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

        // two matrices for each node less the root
        BufferIndexHelper matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
        BufferIndexHelper scaleBufferHelper = new BufferIndexHelper(scaleBufferCount, 0);
        
        //null implies no restrictions
        int[] resourceList = null;
        // TODO
        long preferenceFlags = 0; 
        // TODO
        long requirementFlags = 0;
        
        Beagle beagle = BeagleFactory.loadBeagleInstance(
        tipCount,
        partialBufferHelper.getBufferCount(),
        compactPartialsCount,
        stateCount,
        patternCount,
        eigenBufferHelper.getBufferCount(), // eigenBufferCount
        matrixBufferHelper.getBufferCount() + branchSubstitutionModel.getExtraBufferCount(treeModel),
        categoryCount,
        scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
        resourceList,
        preferenceFlags,
        requirementFlags
);

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
//		PatternList patternList = null;
//		TipStatesModel tipStatesModel = null;
//		boolean useAmbiguities = false;
//		PartialsRescalingScheme rescalingScheme = null;
//		
//		BeagleTreeLikelihood treeLikelihood = new BeagleTreeLikelihood(
//				patternList, //
//				treeModel, //
//				branchSubstitutionModel, //
//				siteRateModel, //
//				branchRateModel, //
//				tipStatesModel, //
//				useAmbiguities, //
//				rescalingScheme //
//		);
		
	} // END: Constructor
	
	 public static void main(String [] args) {
	    	
	    	try {
	    	
	    		int sequenceLength = 10;

	    		// create tree
	    		NewickImporter importer = new NewickImporter("(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
	    		Tree tree =  importer.importTree(null);
	    		TreeModel treeModel = new TreeModel(tree);
	    		
	    		// create site model
	    		GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
	    		
	    		// create Frequency Model
	            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
	            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
	            List<FrequencyModel> frequencyModelList = new ArrayList<FrequencyModel>();
	            frequencyModelList.add(freqModel);
	            
	            // create Epoch Model
	    		Parameter kappa1 = new Parameter.Default(1, 1);
	    		Parameter kappa2 = new Parameter.Default(1, 10);
	            HKY hky1 = new HKY(kappa1, freqModel);
	            HKY hky2 = new HKY(kappa2, freqModel);
	        	List<SubstitutionModel> substModelList = new ArrayList<SubstitutionModel>();
				substModelList.add(hky1);
				substModelList.add(hky2);
				
				Parameter epochTimes = new Parameter.Default(1, 20);
				EpochBranchSubstitutionModel substitutionModel = new EpochBranchSubstitutionModel(
						substModelList, //
						frequencyModelList, //
						epochTimes //
				);
	        	
	    		// create branch rate model
	    		BranchRateModel branchRateModel = new DefaultBranchRateModel();

//	    		 feed to sequence simulator and generate leaves
	    		BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(treeModel, //
	    		substitutionModel, //		
	    		siteRateModel, //
	    		branchRateModel, //
	    		freqModel, //
	    		sequenceLength //
	    		);

//	    		System.out.println(treeSimulator.simulate().toString());

			} catch (Exception e) {
				e.printStackTrace();
			}//END: try-catch block
			
		} // END: main
	
} //END: class
