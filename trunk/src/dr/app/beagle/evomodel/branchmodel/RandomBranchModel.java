package dr.app.beagle.evomodel.branchmodel;

import java.util.LinkedHashMap;
import java.util.List;

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class RandomBranchModel extends AbstractModel implements BranchModel{

    public static final String RANDOM_BRANCH_MODEL = "randomBranchModel";
    private final TreeModel treeModel;
    private final List<SubstitutionModel> substitutionModels;
    
//    private int[] order;
    private LinkedHashMap<NodeRef, Integer> branchAssignmentMap;
    
	public RandomBranchModel(TreeModel treeModel,
            List<SubstitutionModel> substitutionModels) {
		
		super(RANDOM_BRANCH_MODEL);
		
		
		this.treeModel = treeModel;
		this.substitutionModels = substitutionModels;
		
		int nodeCount = treeModel.getNodeCount();
		int nModels = substitutionModels.size();
		
		// randomly decide order, once and for all
		  branchAssignmentMap = new LinkedHashMap<NodeRef, Integer>();
		for (int i = 0; i < nodeCount; i++) {

			NodeRef node = treeModel.getNode(i);
			int branchClass = MathUtils.nextInt(nModels);
			branchAssignmentMap.put(node, branchClass);

		}// END: nodes loop		
		
	}//END: Constructor
	
	@Override
	public Mapping getBranchModelMapping(NodeRef branch) {
		
		final int branchClass = branchAssignmentMap.get(branch);
		
        return new Mapping() {
            public int[] getOrder() {
                return new int[] { branchClass };
            }

            public double[] getWeights() {
                return new double[] { 1.0 };
            }
        };
	}

	@Override
	public List<SubstitutionModel> getSubstitutionModels() {
		return substitutionModels;
	}

	@Override
	public SubstitutionModel getRootSubstitutionModel() {
		int rootClass = branchAssignmentMap.get(treeModel.getRoot());
		return substitutionModels.get(rootClass);
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
		fireModelChanged();
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


}//END: class
