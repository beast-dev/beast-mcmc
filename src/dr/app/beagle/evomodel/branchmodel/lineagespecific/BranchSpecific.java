package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.BranchModel.Mapping;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

public class BranchSpecific extends AbstractBranchRateModel implements BranchModel {

	private Parameter rateParameter;
	
	private Map<NodeRef, Mapping> nodeMap = new HashMap<NodeRef, Mapping>();
	private SubstitutionModel codonModel;
	private boolean setupMapping = true;
	
	private CountableBranchCategoryProvider rateCategories; 
	private TreeModel treeModel;
	
	public BranchSpecific(TreeModel treeModel, SubstitutionModel codonModel,CountableBranchCategoryProvider rateCategories) {
		
		super("");
		
		this.treeModel = treeModel;
		this.codonModel = codonModel;

		this.rateCategories = rateCategories;
		
	}
	
	@Override
	public double getBranchRate(Tree tree, NodeRef node) {
		
		//form product? or draw from binomial for alpha and beta?
		
		
		int rateCategory = rateCategories.getBranchCategory(treeModel, node);
        double value = rateParameter.getParameterValue(rateCategory);
		
		
		return value;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
	}

	@Override
	protected void storeState() {
	}

	@Override
	protected void restoreState() {
	}

	@Override
	protected void acceptState() {
	}

	@Override
	public Mapping getBranchModelMapping(NodeRef branch) {
		
		if(setupMapping ) {
		setupNodeMap(branch);
		}
		
		return nodeMap.get(branch);
	}//END: getBranchModelMapping
	
	
	public void setupNodeMap(NodeRef branch) {
		
		
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
		list.add(codonModel);
		
		return list;
	}

	@Override
	public SubstitutionModel getRootSubstitutionModel() {
		return codonModel;
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

	
}// END: class
