package dr.app.beauti.components.dnds;

import java.util.ArrayList;
import java.util.List;

import dr.app.beauti.options.*;

/**
 * @author Filip Bielejec
 * @version $Id$
 */

public class DnDsComponentOptions implements ComponentOptions {

	final private BeautiOptions options;
	static public final String CODON_PARTITIONED_ROBUST_COUNTING = "codon.partitioned.robust.counting";
	private List<PartitionSubstitutionModel> partitionList = new ArrayList<PartitionSubstitutionModel>();

	public DnDsComponentOptions(final BeautiOptions options) {
		this.options = options;
	}

	public void createParameters(ModelOptions modelOptions) {
		// Do nothing; this is only called at launch
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

	public boolean addPartition(PartitionSubstitutionModel partitionModel) {
		if (!partitionList.contains(partitionModel)) {
			partitionList.add(partitionModel);
		}
		return true; // No error
	}

	public void removePartition(PartitionSubstitutionModel partitionModel) {
		if (partitionList.contains(partitionModel)) {
			partitionList.remove(partitionModel);
		}
	}

	public List<PartitionSubstitutionModel> getPartitionList() {
		return partitionList;
	}

	public BeautiOptions getOptions() {
		return options;
	}

    public boolean doRobustCounting(PartitionSubstitutionModel partitionModel) {
        return partitionList.contains(partitionModel);
    }

}
