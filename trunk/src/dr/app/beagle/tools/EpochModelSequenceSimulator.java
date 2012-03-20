
package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;


import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class EpochModelSequenceSimulator {
	
	private int sequenceLength;
	private Tree tree;
	private SiteRateModel siteModel;
	private BranchRateModel branchRateModel;
	private FrequencyModel freqModel;
	private EpochBranchSubstitutionModel substitutionModel;
	private int categoryCount;
	private int stateCount;

	private static boolean hasAncestralSequence = false;
	private Sequence ancestralSequence;
    
    private double[][] transitionProbabilities;

	EpochModelSequenceSimulator(Tree tree, //
			SiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			EpochBranchSubstitutionModel substitutionModel, //
			int sequenceLength //
	) {

		this.tree = tree;

		this.siteModel = siteModel;
		this.branchRateModel = branchRateModel;
		this.freqModel = freqModel;
		this.substitutionModel = substitutionModel;

		this.sequenceLength = sequenceLength;
		this.stateCount = freqModel.getDataType().getStateCount();
		this.categoryCount = siteModel.getCategoryCount();
		this.transitionProbabilities = new double[categoryCount][stateCount * stateCount];
	} // END: Constructor

	Sequence intArray2Sequence(int[] seq, NodeRef node) {
		StringBuilder sSeq = new StringBuilder();
		for (int i = 0; i < sequenceLength; i++) {
			sSeq.append(freqModel.getDataType().getCode(seq[i]));
		}

		return new Sequence(tree.getNodeTaxon(node), sSeq.toString());
	} // END: intArray2Sequence

	int[] sequence2intArray(Sequence seq) {

		if(seq.getLength() != sequenceLength) {

			throw new RuntimeException("Ancestral sequence length has "
					+ seq.getLength() + " characters " + "expecting "
					+ sequenceLength + " characters");

		}
		
		int array[] = new int[sequenceLength];
		for (int i = 0; i < sequenceLength; i++) {
			array[i] = freqModel.getDataType().getState(
					seq.getChar(i));
		}
		return array;
	}//END: sequence2intArray
	
	void setAncestralSequence(Sequence seq) {
		ancestralSequence = seq;
		hasAncestralSequence = true;
	} // END: setAncestralSequence
	
	/**
	 * Perform the actual sequence generation
	 * @return: Alignment containing randomly generated sequences for the nodes in the
	 * leaves of the tree
	 */
	public Alignment simulate() {
		
		NodeRef root = tree.getRoot();

		double[] categoryProbs = siteModel.getCategoryProportions();
		int[] category = new int[sequenceLength];
		for (int i = 0; i < sequenceLength; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

		int[] seq = new int[sequenceLength];

		if (hasAncestralSequence) {

			seq = sequence2intArray(ancestralSequence);
			
		} else {

			for (int i = 0; i < sequenceLength; i++) {
				seq[i] = MathUtils.randomChoicePDF(freqModel
						.getFrequencies());
			}

		}

		SimpleAlignment alignment = new SimpleAlignment();
		alignment.setReportCountStatistics(false);
		alignment.setDataType(freqModel.getDataType());

		traverse(root, seq, category, alignment);

		return alignment;
	} // END: simulate

	void traverse(NodeRef node, int [] parentSequence, int [] category, SimpleAlignment alignment) {
		
		for (int iChild = 0; iChild < tree.getChildCount(node); iChild++) {
			
			NodeRef child = tree.getChild(node, iChild);
			
			for (int i = 0; i < categoryCount; i++) {
				
				substitutionModel.getTransitionProbabilities(tree, //
						child, //
						i, //
						branchRateModel, //
						siteModel, //
						freqModel, //
						transitionProbabilities[i] //
						);
			
//				EpochBranchSubstitutionModel.printArray(transitionProbabilities[i], transitionProbabilities.length);
				
			}// END: category count loop

        	int [] seq = new int[sequenceLength];
    		double [] cProb = new double[stateCount];
        	for (int i  = 0; i < sequenceLength; i++) {
        		
        		System.arraycopy(transitionProbabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
        		
//        		EpochBranchSubstitutionModel.printArray(cProb, cProb.length);
        		
            	seq[i] = MathUtils.randomChoicePDF(cProb);
            	
        	}

            if (tree.getChildCount(child) == 0) {
            	alignment.addSequence(intArray2Sequence(seq, child));
            }
            
			traverse(tree.getChild(node, iChild), seq, category, alignment);
		}

	} // END: traverse

  public static final String BEAGLE_SEQUENCE_SIMULATOR = "beagleSequenceSimulator";
  public static final String SITE_MODEL ="siteModel";
  public static final String TREE = "tree";
  public static final String EPOCH_BRANCH_SUBSTITUTION_MODEL = "epochBranchSubstitutionModel";
  public static final String REPLICATIONS = "replications";

  public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

      public String getParserName() {
          return BEAGLE_SEQUENCE_SIMULATOR;
      }

      public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		  int replications = xo.getIntegerAttribute(REPLICATIONS);

          Tree tree = (Tree) xo.getChild(Tree.class);
          GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
          BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);
          EpochBranchSubstitutionModel substitutionModel = (EpochBranchSubstitutionModel)xo.getChild(EpochBranchSubstitutionModel.class);
          
          FrequencyModel freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
          
          Sequence ancestralSequence = (Sequence)xo.getChild(Sequence.class);
          
          if (branchRateModel == null) {
        	  branchRateModel = new DefaultBranchRateModel();
          }
          
          EpochModelSequenceSimulator s = new EpochModelSequenceSimulator(tree, //
        		  siteModel, //
        		  branchRateModel, //
        		  freqModel, //
        		  substitutionModel, //
        		  replications //
        		  );

          if(ancestralSequence != null) {
              s.setAncestralSequence(ancestralSequence);
          }
          
          return s.simulate();
      }

      //************************************************************************
      // AbstractXMLObjectParser implementation
      //************************************************************************

      public String getParserDescription() {
          return "An EpochSequenceSimulator that generates random sequences for a given tree, sitemodel and epoch substitution model";
      }

      public Class<Alignment> getReturnType() {
          return Alignment.class;
      }

      public XMLSyntaxRule[] getSyntaxRules() {
          return rules;
      }

      private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
              new ElementRule(Tree.class),
              new ElementRule(GammaSiteRateModel.class),
              new ElementRule(BranchRateModel.class, true),
              new ElementRule(EpochBranchSubstitutionModel.class, true),
              new ElementRule(FrequencyModel.class, true),
              new ElementRule(Sequence.class, true),
              AttributeRule.newIntegerRule(REPLICATIONS)
      };
  };
	
    public static void main(String [] args) {
    	
    	try {
    	
    		int sequenceLength = 10;

    		// create tree
    		NewickImporter importer = new NewickImporter("(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
    		Tree tree =  importer.importTree(null);
    		
    		// create site model
    		GammaSiteRateModel siteModel = new GammaSiteRateModel("siteModel");
    		
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

    		// feed to sequence simulator and generate leaves
    		EpochModelSequenceSimulator treeSimulator = new EpochModelSequenceSimulator(tree, //
    		siteModel, //
    		branchRateModel, //
    		freqModel, //
    		substitutionModel, //
    		sequenceLength //
    		);

    		System.out.println(treeSimulator.simulate().toString());

		} catch (Exception e) {
			e.printStackTrace();
		}//END: try-catch block
		
	} // END: main

} //END: class 
