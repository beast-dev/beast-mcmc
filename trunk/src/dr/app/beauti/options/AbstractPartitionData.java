package dr.app.beauti.options;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.util.TaxonList;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public abstract class AbstractPartitionData {
    protected String fileName;
    protected String name;
    protected TraitData trait;

    protected BeautiOptions options;
    protected PartitionSubstitutionModel model;
    protected PartitionClockModel clockModel;
    protected PartitionTreeModel treeModel;

    protected double meanDistance;

    protected DistanceMatrix distances;

    protected void calculateMeanDistance(Patterns patterns) {
        if (patterns != null) {
            distances = new JukesCantorDistanceMatrix(patterns);
            meanDistance = distances.getMeanDistance();
        } else {
            distances = null;
            meanDistance = 0.0;
        }
    }

    public String getFileName() {
        return fileName;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public TraitData getTrait() {
        return trait;
    }

    public double getMeanDistance() {
        return meanDistance;
    }

    public DistanceMatrix getDistances() {
        return distances;
    }

    public void setPartitionSubstitutionModel(PartitionSubstitutionModel model) {
        this.model = model;
    }

    public PartitionSubstitutionModel getPartitionSubstitutionModel() {
        return this.model;
    }

    public void setPartitionClockModel(PartitionClockModel clockModel) {
        this.clockModel = clockModel;
    }

    public PartitionClockModel getPartitionClockModel() {
        return clockModel;
    }

    public PartitionTreeModel getPartitionTreeModel() {
        return treeModel;
    }

    public void setPartitionTreeModel(PartitionTreeModel treeModel) {
        this.treeModel = treeModel;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.dataPartitions != null && options.dataPartitions.size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public int getTaxonCount() {
        if (getTaxonList() != null) {
            return getTaxonList().getTaxonCount();
        } else {
            // is a trait
            return -1;
        }
    }

    public abstract TaxonList getTaxonList();

    public abstract int getSiteCount();

    public abstract DataType getDataType();

    public abstract String getDataDescription();
}
