package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.options.*;
import dr.app.beauti.types.*;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AncestralStatesComponentOptions implements ComponentOptions {

    public AncestralStatesComponentOptions() {
    }

    public void createParameters(final ModelOptions modelOptions) {
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
    }

    public boolean reconstructAtNodes(final AbstractPartitionData partition) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
        return options.reconstructAtNodes;
    }

    public void setReconstructAtNodes(final AbstractPartitionData partition, boolean reconstructAtNodes) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
         options.reconstructAtNodes = reconstructAtNodes;
    }

    public boolean reconstructAtMRCA(final AbstractPartitionData partition) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
        return options.reconstructAtMRCA;
    }

    public String getMRCAName(final AbstractPartitionData partition) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
        return options.mrcaName;
    }

    public boolean robustCounting(final AbstractPartitionData partition) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
        return options.robustCounting;
    }

    public boolean dNdSRobustCounting(final AbstractPartitionData partition) {
        AncestralStateOptions options = ancestralStateOptionsMap.get(partition);
        return options.dNdSRobustCounting;
    }

    class AncestralStateOptions {
        boolean reconstructAtNodes = false;
        boolean reconstructAtMRCA = false;
        String mrcaName = null;
        boolean robustCounting = false;
        boolean dNdSRobustCounting = false;
    };

    private final Map<AbstractPartitionData, AncestralStateOptions> ancestralStateOptionsMap = new HashMap<AbstractPartitionData, AncestralStateOptions>();

}