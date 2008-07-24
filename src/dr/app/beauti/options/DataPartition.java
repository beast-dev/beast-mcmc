package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DataPartition {

    public DataPartition(String name, String fileName, Alignment alignment) {
        this.name = name;
        this.fileName = fileName;
        this.alignment = alignment;
        this.coding = false;

        Patterns patterns = new Patterns(alignment);
        DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
        meanDistance = distances.getMeanDistance();

    }

    public double getMeanDistance() {
        return meanDistance;
    }

    public String getFileName() {
        return fileName;
    }

    public Alignment getAlignment() {
        return alignment;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setPartitionModel(PartitionModel model) {
        this.model = model;
    }

    public PartitionModel getPartitionModel() {
        return model;
    }

    public boolean isCoding() {
        return coding;
    }

    public void setCoding(boolean coding) {
        this.coding = coding;
    }

    public String toString() {
        return getName();
    }

    private final String fileName;
    private final Alignment alignment;
    private final double meanDistance;

    private String name;
    private boolean coding;


    private PartitionModel model;
}
