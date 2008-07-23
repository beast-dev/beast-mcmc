package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;

/**
 * @author Andrew Rambaut
 */
public class DataPartition {

    public DataPartition(String name, String fileName, Alignment alignment) {
        this.name = name;
        this.fileName = fileName;
        this.alignment = alignment;
        this.coding = false;
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

    private final String fileName;
    private final Alignment alignment;

    private String name;
    private boolean coding;

    private PartitionModel model;
}
