package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;

import beagle.Beagle;
import beagle.BeagleFactory;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
//import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Class for performing random sequence generation for a given site model.
 * Sequences for the leave nodes in the tree are returned as an alignment.
 */

public class BeagleSequenceSimulator {

private TreeModel treeModel;
private int nReplications;
private GammaSiteRateModel siteModel;
private FrequencyModel freqModel;

	/**
	 * Constructor that invokes Beagle
	 */
BeagleSequenceSimulator(TreeModel treeModel, BranchSubstitutionModel branchSubstitutionModel, int nReplications, FrequencyModel freqModel, GammaSiteRateModel siteModel) {
		
		this.treeModel = treeModel;
		this.nReplications = nReplications;
		this.siteModel = siteModel;
		this.freqModel = freqModel;

		int tipCount = treeModel.getExternalNodeCount();
		int nodeCount = treeModel.getNodeCount();
		int eigenCount = branchSubstitutionModel.getEigenCount();
		int internalNodeCount = treeModel.getInternalNodeCount();
		int scaleBufferCount = internalNodeCount + 1;
		
		// TODO: useAmbiguities ?
		int compactPartialsCount = tipCount;
		
		// TODO: do we need this when simulating?
		int patternCount = nReplications;
		
		// TODO: don't really need a siteModel to get the stateCount
		// too bad the siteModel only functions as a container with respect to the substitutionModel
		// parse from BranchSubstitutionModel ?
		int stateCount = freqModel.getDataType().getStateCount();
		
		// TODO: don't like this categoryCount though, does require a siteModel
		int categoryCount = siteModel.getCategoryCount();
		
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

		
	}//END: Constructor
	
	public Alignment simulate() {

		NodeRef root = treeModel.getRoot();

		double[] categoryProbs = siteModel.getCategoryProportions();
		int[] category = new int[nReplications];
		for (int i = 0; i < nReplications; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

		int[] seq = new int[nReplications];

		for (int i = 0; i < nReplications; i++) {
			seq[i] = MathUtils.randomChoicePDF(freqModel.getFrequencies());
		}

		SimpleAlignment alignment = new SimpleAlignment();
		alignment.setReportCountStatistics(false);
		alignment.setDataType(freqModel.getDataType());

//		traverse(root, seq, category, alignment);

		return alignment;
	} // END: simulate

	public static void main(String[] args) {

		try {

			
		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	} // END: main

//	public static final String BEAGLE_SEQUENCE_SIMULATOR = "beagleSequenceSimulator";
//    public static final String SITE_MODEL = SiteModel.SITE_MODEL;
//    public static final String TREE = "tree";
//    public static final String EPOCH_BRANCH_SUBSTITUTION_MODEL = "epochBranchSubstitutionModel";
//    public static final String REPLICATIONS = "replications";

//    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
//
//        public String getParserName() {
//            return EPOCH_SEQUENCE_SIMULATOR;
//        }
//
//        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//
//        	int nReplications = xo.getIntegerAttribute(REPLICATIONS);
//
//            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
//            GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
//            EpochBranchSubstitutionModel epochModel = (EpochBranchSubstitutionModel)xo.getChild(EpochBranchSubstitutionModel.class);
//            
//            Sequence ancestralSequence = (Sequence)xo.getChild(Sequence.class);
//            
//            EpochBranchSubstitutionModelSimulator s = new EpochBranchSubstitutionModelSimulator(treeModel, epochModel, siteModel, nReplications);
//
//            if(ancestralSequence != null) {
//                s.setAncestralSequence(ancestralSequence);
//            }
//            
//            return s.simulate();
//        }
//
//        //************************************************************************
//        // AbstractXMLObjectParser implementation
//        //************************************************************************
//
//        public String getParserDescription() {
//            return "An EpochSequenceSimulator that generates random sequences for a given tree, sitemodel and epoch substitution model";
//        }
//
//        public Class<Alignment> getReturnType() {
//            return Alignment.class;
//        }
//
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return rules;
//        }
//
//        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//                new ElementRule(TreeModel.class),
//                new ElementRule(GammaSiteRateModel.class),
//                new ElementRule(BranchRateModel.class, true),
//                new ElementRule(EpochBranchSubstitutionModel.class, true),
//                new ElementRule(Sequence.class, true),
//                AttributeRule.newIntegerRule(REPLICATIONS)
//        };
//    };

//	/** generate simple site model, for testing purposes **/
//	static GammaSiteRateModel getDefaultGammaSiteRateModel() {
//		
//		List<FrequencyModel> frequencyModelList = new ArrayList<FrequencyModel>();
//		List<SubstitutionModel> substModelList = new ArrayList<SubstitutionModel>();
//		Parameter epochTransitionTimes = new Parameter.Default(1, 20);
//		
//		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
//				0.25, 0.25 });
//		FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
//		
//		Parameter kappa1 = new Parameter.Default(1, 1);
//		Parameter kappa2 = new Parameter.Default(10, 1);
//		
//		HKY hky1 = new HKY(kappa1, freqModel);
//		HKY hky2 = new HKY(kappa2, freqModel);
//		
//		substModelList.add((SubstitutionModel) hky1);
//		substModelList.add((SubstitutionModel) hky2);
//		frequencyModelList.add(freqModel);
//		
//		GammaSiteRateModel gsrm = new GammaSiteRateModel("dupa");
//		gsrm.setSubstitutionModel(hky1);
//
//		return gsrm;
//	} // END: getDefaultGammaSiteRateModel
	
} // class SequenceSimulator
