package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DataPartition {

    public DataPartition(String name, String fileName, Alignment alignment) {
        this(name, fileName, alignment, -1, -1);
    }

    public DataPartition(String name, String fileName, Alignment alignment, int fromSite, int toSite) {
        this.name = name;
        this.fileName = fileName;
        this.alignment = alignment;
        this.coding = false;

        this.fromSite = fromSite;
        this.toSite = toSite;

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

    public int getSiteCount() {
        return getFromSite() - getToSite() + 1;
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

    private PartitionModel model;

    public void addPatternListId(String s) {
        patternListIds.add(s);

    }

    public String getPatternListId(int i) {
        return patternListIds.get(i);
    }

    List<String> patternListIds = new ArrayList<String>();

}
