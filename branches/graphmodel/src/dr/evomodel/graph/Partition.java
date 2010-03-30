package dr.evomodel.graph;

import dr.evolution.alignment.SiteList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/*
 * Refers to a range of sites in a SiteList
 * Used for mapping lineage of sites in a graph model
 */
public class Partition extends AbstractModel {
	public static final  String PARTITION = "partition";
	int number = -1;

	SiteList siteList = null;
	int leftSite = -1;
	int rightSite = -1;	
	SiteModel siteModel = null;
	BranchRateModel branchRateModel = null;
	
	public Partition(){
		super(PARTITION);
	}
	public Partition(SiteList siteList){
		this(siteList, 0, siteList.getSiteCount());
	}
	public Partition(SiteList siteList, int leftSite, int rightSite){
		super(PARTITION);
		this.siteList = siteList;
		this.leftSite = leftSite;
		this.rightSite = rightSite;
	}
	
	public SiteList getSiteList() {
		return siteList;
	}
	public void setSiteList(SiteList siteList) {
		this.siteList = siteList;
	}
	public int getLeftSite() {
		return leftSite;
	}
	public void setLeftSite(int leftSite) {
		this.leftSite = leftSite;
	}
	public int getRightSite() {
		return rightSite;
	}
	public void setRightSite(int rightSite) {
		this.rightSite = rightSite;
	}
	
	public int getNumber() {
	    return number;
	}
	
	public void setNumber(int n) {
	    number = n;
	}

	public SiteModel getSiteModel() {
		return siteModel;
	}
	public void setSiteModel(SiteModel siteModel) {
		if(this.siteModel != null)
			this.siteModel.removeModelListener(this);
		this.siteModel = siteModel;
		if(this.siteModel != null)
			siteModel.addModelListener(this);
	}
	public BranchRateModel getBranchRateModel() {
		return branchRateModel;
	}
	public void setBranchRateModel(BranchRateModel branchRateModel) {
		if(this.branchRateModel != null)
			this.branchRateModel.removeModelListener(this);
		this.branchRateModel = branchRateModel;
		if(this.branchRateModel != null)
			branchRateModel.addModelListener(this);
	}

	public void copyPartition(Partition p){
	       setLeftSite(p.getLeftSite());
	       setRightSite(p.getRightSite());
	       setSiteList(p.getSiteList());
	       setSiteModel(p.getSiteModel());
	       setBranchRateModel(p.getBranchRateModel());
	}

	SiteList storedSiteList = null;
	int storedLeftSite = -1;
	int storedRightSite = -1;	
	SiteModel storedSiteModel = null;
	BranchRateModel storedBranchRateModel = null;

	//
	// BEGIN AbstractModel INTERFACE
	//
	
	@Override
	protected void acceptState() {
	}

	@Override
	protected void restoreState() {
		setSiteList(storedSiteList);
		setLeftSite(storedLeftSite);
		setRightSite(storedRightSite);	
		setSiteModel(storedSiteModel);
		setBranchRateModel(storedBranchRateModel);
	}

	@Override
	protected void storeState() {
		storedSiteList = siteList;
		storedLeftSite = leftSite;
		storedRightSite = rightSite;
		storedSiteModel = siteModel;
		storedBranchRateModel = branchRateModel;
	}
	
	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// forward the change in the submodel (branchRateModel or siteModel)
		fireModelChanged(object);
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// this shouldn't happen
	}
	//
	// END AbstractModel INTERFACE
	//
}
