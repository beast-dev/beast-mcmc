package dr.app.beagle.evomodel.branchmodel;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math.random.MersenneTwister;

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.GY94CodonModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.bss.Utils;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Parameter.Default;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */

@SuppressWarnings("serial")
public class RandomBranchModel extends AbstractModel implements BranchModel {

	public static final String RANDOM_BRANCH_MODEL = "randomBranchModel";
	private final TreeModel treeModel;
	private GY94CodonModel baseSubstitutionModel;

	private LinkedList<SubstitutionModel> substitutionModels;

	// private int[] order;
	private LinkedHashMap<NodeRef, Integer> branchAssignmentMap;

	private static MersenneTwister random = new MersenneTwister(
			MathUtils.nextLong());

	public RandomBranchModel(TreeModel treeModel,
			GY94CodonModel baseSubstitutionModel) {

		super(RANDOM_BRANCH_MODEL);

		this.treeModel = treeModel;
		this.baseSubstitutionModel = baseSubstitutionModel;

		setup();

	}// END: Constructor

	private void setup() {

		// TODO: parse
		double stdev = 1.0;
		double mean = 0.0;

		DataType dataType = baseSubstitutionModel.getDataType();
		FrequencyModel freqModel = baseSubstitutionModel.getFrequencyModel();
		Parameter kappaParameter = new Parameter.Default("kappa", 1, baseSubstitutionModel.getKappa());
		
		substitutionModels = new LinkedList<SubstitutionModel>();
		branchAssignmentMap = new LinkedHashMap<NodeRef, Integer>();

		int branchClass = 0;
		for (NodeRef node : treeModel.getNodes()) {

			double time = treeModel.getNodeHeight(node);
			double baseOmega = baseSubstitutionModel.getOmega();
			double epsilon = (random.nextGaussian() * stdev + mean);

			double value = baseOmega * time + epsilon;
			Parameter omegaParameter = new Parameter.Default("omega", 1, value);
			GY94CodonModel gy94 = new GY94CodonModel((Codons) dataType, omegaParameter, kappaParameter, freqModel);

			substitutionModels.add(gy94);
			branchAssignmentMap.put(node, branchClass);
			branchClass++;
		}// END: nodes loop

	}// END: setup

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

}// END: class
