package dr.app.bss;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> implements Serializable {

	public MutableTaxonList taxonList = new Taxa();
	public int simulationsCount = 1;
	
	public LinkedList<File> treeFileList = new LinkedList<File>();
	public LinkedList<Integer> taxaCounts = new LinkedList<Integer>();

	// do not serialize this two
	public transient boolean setSeed = false;
	public transient long startingSeed;
	
	public PartitionDataList() {
		super();
		startingSeed = System.currentTimeMillis();
	}// END: Constructor

}// END:class
