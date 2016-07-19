package dr.evomodel.branchmodel;

import dr.evomodelxml.branchmodel.BranchAssignmentModelParser;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

@SuppressWarnings("serial")
public class BranchAssignmentModel extends AbstractModel implements BranchModel {

	public static final String BRANCH_ASSIGNMENT_MODEL = "branchAssignmentModel";
	
	private TreeModel treeModel;
	private final String annotation;
	private final LinkedHashMap<Integer, SubstitutionModel> modelIndexMap;
	private final SubstitutionModel baseModel;

	private Integer baseModelIndex;
	private LinkedHashMap<NodeRef, Integer> branchAssignmentMap;
	private LinkedList<SubstitutionModel> substitutionModels;

	public BranchAssignmentModel(
			TreeModel treeModel, //
			String annotation, //
			LinkedHashMap<Integer, SubstitutionModel> modelIndexMap, //
			SubstitutionModel baseModel//
	) {

		super(BRANCH_ASSIGNMENT_MODEL);

		this.treeModel = treeModel;
		this.annotation = annotation;
		this.modelIndexMap = modelIndexMap;
		this.baseModel = baseModel;

		this.substitutionModels = new LinkedList<SubstitutionModel>();
		this.branchAssignmentMap = new LinkedHashMap<NodeRef, Integer>();
		// base model comes last
		this.baseModelIndex = modelIndexMap.size();
		
		setup();

	}// END: Constructor

	private void setup() {
		
		// for (int i = 0; i < modelIndexMap.size() + 1; i++) {
		// substitutionModels.add(null);
		// }
		
//	try {	
//		
//		File file = new File("/home/filip/Dropbox/BeagleSequenceSimulator/branchSpecificSimulations/annotated_tree.nexus");
//		BufferedReader reader;
//		
//		reader = new BufferedReader(new FileReader(file));
//		NexusImporter importer = new NexusImporter(reader);
//		Tree tree = importer.importTree(null);
//		this.treeModel = new TreeModel(tree);
//		
//	} catch ( Exception e) {
//		e.printStackTrace();
//	} 
		
		for (NodeRef node : this.treeModel.getNodes()) {
			if (!treeModel.isRoot(node)) {

				Integer modelIndex = Integer.MAX_VALUE;
				SubstitutionModel model = null;
				Object nodeAttribute = treeModel.getNodeAttribute(node,
						annotation);

				if (nodeAttribute == null) {

					System.out
							.println("Attribute "
									+ annotation
									+ " missing from node. Using base model as branch model.");

					modelIndex = this.baseModelIndex;
					model = this.baseModel;

				} else {

					modelIndex = (Integer) nodeAttribute;
					model = this.modelIndexMap.get(modelIndex);

				}

				branchAssignmentMap.put(node, modelIndex);

				// if (substitutionModels.get(modelIndex) == null) {
				// substitutionModels.set(modelIndex, model);
				// }

				substitutionModels.add(model);

			}// END: root check
		}// END: nodes loop
		
	}//END: setup
	
	@Override
	public Mapping getBranchModelMapping(NodeRef branch) {

		final int modelIndex = branchAssignmentMap.get(branch);

		return new Mapping() {
			public int[] getOrder() {
				return new int[] { modelIndex };
			}

			public double[] getWeights() {
				return new double[] { 1.0 };
			}
		};
	}// END: getBranchModelMapping

	@Override
	public List<SubstitutionModel> getSubstitutionModels() {
		return substitutionModels;
	}// END: getSubstitutionModels

	@Override
	public SubstitutionModel getRootSubstitutionModel() {

		Object nodeAttribute = treeModel.getNodeAttribute(treeModel.getRoot(),
				BranchAssignmentModelParser.ANNOTATION_VALUE);
		SubstitutionModel model = null;

		if (nodeAttribute == null) {

			model = this.baseModel;

		} else {

			Integer modelIndex = (Integer) nodeAttribute;
			model = this.modelIndexMap.get(modelIndex);

		}

		return model;
	}// END: getRootSubstitutionModel

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
	protected void handleVariableChangedEvent(@SuppressWarnings("rawtypes") Variable variable, int index,
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

}// END: class
