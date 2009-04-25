package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DataPartition {

    public DataPartition(String name, String fileName, Alignment alignment) {
        this(name, fileName, alignment, -1, -1, 1);
    }

    public DataPartition(String name, String fileName, Alignment alignment, int fromSite, int toSite, int every) {
        this.name = name;
        this.fileName = fileName;
        this.alignment = alignment;
        this.coding = false;

        this.fromSite = fromSite;
        this.toSite = toSite;
        this.every = every;

//        Patterns patterns = new Patterns(alignment);
//        DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
//        meanDistance = distances.getMeanDistance();
        meanDistance = 0.0;

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

    public void setPartitionTree(final PartitionTree tree) {
        this.tree = tree;
    }

    public PartitionTree getPartitionTree() {
        return tree;
    }

    public boolean isCoding() {
        return coding;
    }

    public void setCoding(boolean coding) {
        this.coding = coding;
    }

    public int getFromSite() {
        return fromSite;
    }

    public int getToSite() {
        return toSite;
    }

    public int getEvery() {
        return every;
    }

    public int getSiteCount() {
        int from = getFromSite();
        if (from < 1) {
            from = 1;
        }
        int to = getToSite();
        if (to < 1) {
            to = alignment.getSiteCount();
        }
        return (to - from + 1) / every;
    }

    public String toString() {
        return getName();
    }

    private final String fileName;
    private final Alignment alignment;
    private final double meanDistance;

    private String name;
    private boolean coding;

    private int fromSite;
    private int toSite;
    private int every = 1;


    private PartitionModel model;
    private PartitionTree tree;
}
