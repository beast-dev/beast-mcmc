package dr.app.bss;

import java.io.File;
import java.util.ArrayList;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public MutableTaxonList taxonList = new Taxa();

	public int simulationsCount = 1;

	public File treesFilename = null;

//	public LinkedHashSet<File> treeFileList = new LinkedHashSet<File>();
	public ArrayList<File> treeFileList = new ArrayList<File>();
	
	public PartitionDataList() {
		super();
	}// END: Constructor

}// END:class
