package dr.app.beauti.options;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionPattern extends AbstractPartitionData { // microsatellite

    private final Patterns patterns;

    public PartitionPattern(BeautiOptions options, String name, String fileName, Patterns patterns) {
        this.options = options;
        this.name = name;
        this.fileName = fileName;
        this.patterns = patterns;

        this.trait = null;

        calculateMeanDistance(patterns);
    }

    public Patterns getPatterns() {
        return patterns;
    }

    public TaxonList getTaxonList() {
        return getPatterns();  
    }

    public int getSiteCount() {
        return 1;
    }

    public DataType getDataType() {
        if (patterns != null) {
            return patterns.getDataType();
        } else {
            return trait.getDataType();
        }
    }

    public String getDataDescription() {
        if (patterns != null) {
            return patterns.getDataType().getDescription();
        } else {
            return trait.getTraitType().toString();
        }
    }
}
