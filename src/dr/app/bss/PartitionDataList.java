package dr.app.bss;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;

import dr.evolution.util.Taxa;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> implements Serializable {

	public int simulationsCount = 1;
	public boolean useParallel = false;
	
	//List of all Taxa displayed in Taxa Panel
	public Taxa allTaxa = new Taxa();
	LinkedList<TreesTableRecord> recordsList = new LinkedList<TreesTableRecord>();
	
	// do not serialize this two
	public transient boolean setSeed = false;
	public transient long startingSeed;
	
	public PartitionDataList() {
		super();
		startingSeed = System.currentTimeMillis();
	}// END: Constructor

}// END:class
