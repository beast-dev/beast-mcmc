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
        super(options, name, fileName);
        this.patterns = patterns;

        this.traits = null;

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
            throw new RuntimeException("patterns should not be null");
        }
    }

    public String getDataDescription() {
        if (patterns != null) {
            return patterns.getDataType().getDescription();
        } else {
            throw new RuntimeException("patterns should not be null");
        }
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionPattern().size() > 1) { // getPartitionPattern() already excludes traits and PartitionData
            prefix += getName() + ".";
        }
        return prefix;
    }
}
