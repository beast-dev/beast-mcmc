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

    private List<PartitionSubstitutionModel> partitionList = new ArrayList<PartitionSubstitutionModel>();
    
}
