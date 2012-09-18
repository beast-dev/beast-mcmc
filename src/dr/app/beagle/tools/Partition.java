package dr.app.beagle.tools;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.evolution.sequence.Sequence;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class Partition {

	public int from;
	public int to; 
	public int every;
	
	public TreeModel treeModel;
	public BranchSubstitutionModel branchSubstitutionModel;
	public GammaSiteRateModel siteModel;
	public BranchRateModel branchRateModel;
	public FrequencyModel freqModel;
	
	public boolean hasAncestralSequence = false;
	public Sequence ancestralSequence = null;

//	public static final ThreadLocal threadLocal = new ThreadLocal();
	
	public Partition(
			TreeModel treeModel, //
			BranchSubstitutionModel branchSubstitutionModel,
			GammaSiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int from, //
			int to, //
			int every //
	) {
		
		this.treeModel = treeModel;
		this.siteModel = siteModel;
		this.freqModel = freqModel;
		this.branchSubstitutionModel = branchSubstitutionModel;
		this.branchRateModel = branchRateModel;
		
		this.from = from;
		this.to = to;
		this.every = every;
		
	}//END: Constructor
	
	public void setAncestralSequence(Sequence ancestralSequence) {
		this.ancestralSequence = ancestralSequence;
		this.hasAncestralSequence = true;
	}// END: setAncestralSequence
	
	public int getPartitionSiteCount() {
		int partitionSiteCount = ((to - from) / every) + 1;
		return partitionSiteCount;
	}// END: getPartitionSiteCount
	
}//END: class
