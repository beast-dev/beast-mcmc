package dr.app.beauti.components.dollo;

import dr.app.beauti.options.*;
import dr.app.beauti.types.PriorScaleType;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @version $Id$
 */

public class DolloComponentOptions implements ComponentOptions {

	final private BeautiOptions options;
	static public final String CODON_PARTITIONED_ROBUST_COUNTING = "codon.partitioned.robust.counting";
	private List<PartitionSubstitutionModel> partitionList = new ArrayList<PartitionSubstitutionModel>();

	public DolloComponentOptions(final BeautiOptions options) {
		this.options = options;
	}

	public void createParameters(ModelOptions modelOptions) {
        modelOptions.createNonNegativeParameterUniformPrior("treeModel.tipDates", "date of specified tips",
                PriorScaleType.TIME_SCALE, 1.0, 0.0, Double.MAX_VALUE);

        modelOptions.createScaleOperator("treeModel.tipDates", modelOptions.demoTuning, 3.0);		
	}

	public void selectOperators(ModelOptions modelOptions, List<Operator> ops) {
		// Do nothing
	}

	public void selectParameters(ModelOptions modelOptions,
			List<Parameter> params) {
		// Do nothing
	}

	public void selectStatistics(ModelOptions modelOptions,
			List<Parameter> stats) {
		// Do nothing
	}

	public boolean addPartition(PartitionSubstitutionModel partition) {
		if (!partitionList.contains(partition)) {
			partitionList.add(partition);
		}
		return true; // No error
	}

	public void removePartition(PartitionSubstitutionModel model) {
		if (partitionList.contains(model)) {
			partitionList.remove(model);
		}
	}

	public List<PartitionSubstitutionModel> getPartitionList() {
		return partitionList;
	}

	public BeautiOptions getOptions() {
		return options;
	}

    public boolean doStochasticDollo(PartitionSubstitutionModel model) {
        return false;
    }
}