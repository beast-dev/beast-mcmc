package dr.app.bss;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> implements Serializable {

	public MutableTaxonList taxonList = new Taxa();
	public int simulationsCount = 1;
	public File treesFilename = null;
	public ArrayList<File> treeFileList = new ArrayList<File>();
	public ArrayList<Integer> taxaCounts = new ArrayList<Integer>();
	
	public PartitionDataList() {
		super();
	}// END: Constructor

}// END:class
