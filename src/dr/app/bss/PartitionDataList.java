package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public static final String VERSION = "0.0.1";
	public static final String DATE_STRING = "2012";
	
	public TaxonList taxonList = new Taxa();
	//TODO: last partition should have this number as TO value
	public int sequenceLength = 1000;	
	public int simulationsNumber = 1;
	
	public HashMap<File, TreeModel> forestMap = new HashMap<File, TreeModel>();
	
	public PartitionDataList() {
		super();
	}// END: Constructor
	
//	public int getSequenceLength() {
//		return sequenceLength;
//	}
//	
//	public void setSequenceLength(int sequenceLength) {
//		this.sequenceLength=sequenceLength;
//	}
	
}//END:class
