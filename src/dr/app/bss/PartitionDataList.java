package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public static final String VERSION = "0.0.1";
	public static final String DATE_STRING = "2012";
	
	public MutableTaxonList taxonList = new Taxa();
	
	//TODO: last partition should have this number as TO value
	public int siteCount = 1000;	
	public int simulationsCount = 1;
	
	public LinkedHashMap<File, TreeModel> forestMap = new LinkedHashMap<File, TreeModel>();
	
	public File treesFilename = null;
	
	public PartitionDataList() {
		super();
	}// END: Constructor
	
}// END:class
