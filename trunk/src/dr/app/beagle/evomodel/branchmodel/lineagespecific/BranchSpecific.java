package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.BranchModel.Mapping;
import dr.app.beagle.evomodel.substmodel.AbstractCodonModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.MG94CodonModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

public class BranchSpecific 
extends AbstractModel
//extends AbstractBranchRateModel 
implements BranchModel {

	private Parameter rateParameter;
	
	private Map<NodeRef, Mapping> nodeMap = new HashMap<NodeRef, Mapping>();
	private SubstitutionModel substModel;
	private boolean setupMapping = true;
	
	private CountableBranchCategoryProvider rateCategories; 
	private TreeModel treeModel;
	private FrequencyModel freqModel;
	
	public BranchSpecific(TreeModel treeModel, SubstitutionModel substModel,CountableBranchCategoryProvider rateCategories) {
		
		super("");
		
		this.treeModel = treeModel;
		this.substModel = substModel;
		this.rateCategories = rateCategories;
		
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
		
		int rateCategory = rateCategories.getBranchCategory(treeModel, branch);
        double value = rateParameter.getParameterValue(rateCategory);
        
        Parameter alphaParameter = new Parameter.Default(1,   ((MG94CodonModel) substModel).getAlpha() * value);
        Parameter betaParameter = new Parameter.Default(1,   ((MG94CodonModel) substModel).getBeta() * value);
        
//		codonModel = new MG94CodonModel(codonDataType, alphaParameter, betaParameter, freqModel);
		
        
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
		// TODO Auto-generated method stub
		return null;
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

	
}// END: class
