package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public MutableTaxonList taxonList = new Taxa();

//	public int siteCount = 1000;
	public int simulationsCount = 1;

	public LinkedHashSet<File> forestList = new LinkedHashSet<File>();
	public File treesFilename = null;

	public PartitionDataList() {
		super();
	}// END: Constructor

}// END:class
