package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

//@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	private static final long serialVersionUID = -7703343661531335551L;

	public MutableTaxonList taxonList = new Taxa();

	public int simulationsCount = 1;

	public LinkedHashSet<File> forestList = new LinkedHashSet<File>();
	public File treesFilename = null;

	//TODO: when panel is ready this comes along
	public LinkedHashSet<File> treeFileList = new LinkedHashSet<File>();
	
	public PartitionDataList() {
		super();
	}// END: Constructor

}// END:class
