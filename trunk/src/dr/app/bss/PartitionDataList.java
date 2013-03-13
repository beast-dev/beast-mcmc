package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public MutableTaxonList taxonList = new Taxa();

	// TODO: last partition should have this number as TO value
	public int siteCount = 1000;
	public int simulationsCount = 1;

	public LinkedList<File> forestList = new LinkedList<File>();
	public File treesFilename = null;

	public PartitionDataList() {
		super();
	}// END: Constructor

}// END:class
