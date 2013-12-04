package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.branchmodel.BranchModel.Mapping;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.AbstractCodonModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.MG94CodonModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

public class BranchSpecific 
extends AbstractModel
//extends AbstractBranchRateModel 
implements BranchModel {

	private Parameter uParameter;
	
	private Map<NodeRef, Mapping> nodeMap = new HashMap<NodeRef, Mapping>();
	private SubstitutionModel substModel;
	private boolean setupMapping = true;
	
	private CountableBranchCategoryProvider uCategories; 
	private TreeModel treeModel;
	private FrequencyModel freqModel;
	
	public BranchSpecific(TreeModel treeModel, SubstitutionModel substModel, CountableBranchCategoryProvider uCategories
			, Parameter uParameter
			) {
		
		super("");
		
		this.treeModel = treeModel;
		this.substModel = substModel;
		this.uCategories = uCategories;
		
		
		this.uParameter = uParameter;
		
		this.freqModel = substModel.getFrequencyModel();
		
		
		
		
	}
	
//	public BranchSpecific(TreeModel treeModel, SubstitutionModel codonModel,CountableBranchCategoryProvider rateCategories) {
//		
//		super("");
//		
//		this.treeModel = treeModel;
//		this.codonModel = codonModel;
//
//		this.rateCategories = rateCategories;
//		
//	}
//	
//	@Override
//	public double getBranchRate(Tree tree, NodeRef node) {
//		
//		//form product? or draw from binomial for alpha and beta?
//		
//		int rateCategory = rateCategories.getBranchCategory(treeModel, node);
//        double value = rateParameter.getParameterValue(rateCategory);
//		
//		
//		return value;
//	}
//
//	@Override
//	protected void handleModelChangedEvent(Model model, Object object, int index) {
//	}
//
//	@Override
//	protected void handleVariableChangedEvent(Variable variable, int index,
//			ChangeType type) {
//	}
//
//	@Override
//	protected void storeState() {
//	}
//
//	@Override
//	protected void restoreState() {
//	}
//
//	@Override
//	protected void acceptState() {
//	}

	@Override
	public Mapping getBranchModelMapping(NodeRef branch) {
		
		if(setupMapping ) {
		setupNodeMap(branch);
		}
		
		return nodeMap.get(branch);
	}//END: getBranchModelMapping
	
	
	public void setupNodeMap(NodeRef branch) {
		
		//form product? or draw from binomial for alpha and beta?
		
		int uCategory = uCategories.getBranchCategory(treeModel, branch);
        double uValue = uParameter.getParameterValue(uCategory);
        
        Parameter alphaParameter = new Parameter.Default(1,   ((MG94CodonModel) substModel).getAlpha() * uValue);
        Parameter betaParameter = new Parameter.Default(1,   ((MG94CodonModel) substModel).getBeta() * uValue);
        
        substModel = new MG94CodonModel(Codons.UNIVERSAL, alphaParameter, betaParameter, freqModel);
		
        getSubstitutionModels();
        
		nodeMap.put(branch, new Mapping() {
			@Override
			public int[] getOrder() {
				return new int[]{ 0 };
			}

			@Override
			public double[] getWeights() {
				return new double[] { 1.0 };
			}
		});
	}// END: setupNodeMap

	@Override
	public List<SubstitutionModel> getSubstitutionModels() {
		List<SubstitutionModel> list = new ArrayList<SubstitutionModel>();
		list.add( substModel);
		
		return list;
	}

	@Override
	public SubstitutionModel getRootSubstitutionModel() {
		return substModel;
	}

	@Override
	public FrequencyModel getRootFrequencyModel() {
		return getRootSubstitutionModel().getFrequencyModel();
	}

	@Override
	public boolean requiresMatrixConvolution() {
		return false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void storeState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void restoreState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub
		
	}

	 public static void main(String[] args) {
		 
	try {	
		
		 // the seed of the BEAST
         MathUtils.setSeed(666);
		 
         // create tree
         NewickImporter importer = new NewickImporter(
                 "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
         TreeModel tree = new TreeModel(importer.importTree(null));
		 
      // create site model
         GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                 "siteModel");

         // create branch rate model
         BranchRateModel branchRateModel = new DefaultBranchRateModel();
         
         int sequenceLength = 10;
         ArrayList<Partition> partitionsList = new ArrayList<Partition>();

         // create Frequency Model
       Parameter freqs = new Parameter.Default(new double[]{
     	0.0163936, // 
     	0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
         0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
     	0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
     	0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
     	0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344 //
     	});
         FrequencyModel freqModel = new FrequencyModel(Codons.UNIVERSAL,
                 freqs);

         // create substitution model
         Parameter alpha = new Parameter.Default(1, 10);
         Parameter beta = new Parameter.Default(1, 5);
         MG94CodonModel mg94 = new MG94CodonModel(Codons.UNIVERSAL, alpha, beta, freqModel);
         
         HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(mg94);

         // create partition
         Partition partition1 = new Partition(tree, //
                 substitutionModel,//
                 siteRateModel, //
                 branchRateModel, //
                 freqModel, //
                 0, // from
                 sequenceLength - 1, // to
                 1 // every
         );

		partitionsList.add(partition1);

		// feed to sequence simulator and generate data
		BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
				partitionsList);

		Alignment alignment = simulator.simulate(false);

		ConvertAlignment convert = new ConvertAlignment(Nucleotides.INSTANCE,
				GeneticCode.UNIVERSAL, alignment);
		
		Parameter uParam = new Parameter.Default(1, 0.0001);
		CountableBranchCategoryProvider uCateg = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(tree, uParam);
		
		BranchSpecific branchSpecific = new BranchSpecific(tree, mg94, uCateg, uParam);
		
		BeagleTreeLikelihood nbtl = new BeagleTreeLikelihood(convert, //
				tree, //
				branchSpecific, //
				siteRateModel, //
				branchRateModel, //
				null, //
				false, //
				PartialsRescalingScheme.DEFAULT);

		System.out.println("likelihood = " + nbtl.getLogLikelihood());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}// END: main
	
}// END: class
