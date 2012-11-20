package dr.app.bss;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> {

	public static final String VERSION = "0.0.1";
	public static final String DATE_STRING = "2012";
	
	public TaxonList taxonList = new Taxa();
	public int sequenceLength = 1000;	
	
	public List<File> treeFilesList = new ArrayList<File>();
	public List<TreeModel> treeModelsList= new ArrayList<TreeModel>();
	
//	HashMap<File, TreeModel> forestMap = new HashMap<File, TreeModel>();
	
	public PartitionDataList() {
		super();
//		treeFilesList.add(null);
	}
	
//	public int getSequenceLength() {
//		return sequenceLength;
//	}
//	
//	public void setSequenceLength(int sequenceLength) {
//		this.sequenceLength=sequenceLength;
//	}
	
}//END:class
