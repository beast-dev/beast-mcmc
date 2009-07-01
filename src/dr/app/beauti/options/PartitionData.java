package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionData {

    private final String fileName;
    private final Alignment alignment;
    private final double meanDistance;

    private String name;
    private boolean coding;

    private int fromSite;
    private int toSite;
    private int every = 1;

    private PartitionSubstitutionModel model;
    private PartitionClockModel clockModel;
    private PartitionTreeModel treeModel;
    
//    public List<Tree> userTrees = new ArrayList<Tree>(); // a set of starting tree loaded from NEXUS file
    
	public PartitionData(String name, String fileName, Alignment alignment) {
        this(name, fileName, alignment, -1, -1, 1);
    }

    public PartitionData(String name, String fileName, Alignment alignment, int fromSite, int toSite, int every) {
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

    public void setPartitionSubstitutionModel(PartitionSubstitutionModel model) {
        this.model = model;
    }

    public PartitionSubstitutionModel getPartitionSubstitutionModel() {
        return model;
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
    
    public int getNumOfTaxa() {
    	int n = alignment.getSequenceCount();
    	
    	if (n > 0) {
    		return n;
    	} else {
    		return 0;
    	}
    }

    public String toString() {
        return getName();
    }

}
