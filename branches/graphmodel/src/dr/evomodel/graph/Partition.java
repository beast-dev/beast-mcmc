package dr.evomodel.graph;

import dr.evolution.alignment.SiteList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;

/*
 * Refers to a range of sites in a SiteList
 * Used for mapping lineage of sites in a graph model
 */
public class Partition {
	SiteList siteList = null;
	int leftSite = -1;
	int rightSite = -1;
	int number = -1;
	
	SiteModel siteModel = null;
	BranchRateModel branchRateModel = null;
	FrequencyModel frequencyModel = null;
	
	public Partition(){}
	public Partition(SiteList siteList){
		this(siteList, 0, siteList.getSiteCount());
	}
	public Partition(SiteList siteList, int leftSite, int rightSite){
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
		this.siteModel = siteModel;
	}
	public BranchRateModel getBranchRateModel() {
		return branchRateModel;
	}
	public void setBranchRateModel(BranchRateModel branchRateModel) {
		this.branchRateModel = branchRateModel;
	}
	public FrequencyModel getFrequencyModel() {
		return frequencyModel;
	}
	public void setFrequencyModel(FrequencyModel frequencyModel) {
		this.frequencyModel = frequencyModel;
	}
	public void copyPartition(Partition p){
	       setLeftSite(p.getLeftSite());
	       setRightSite(p.getRightSite());
	       setSiteList(p.getSiteList());
	       setSiteModel(p.getSiteModel());
	       setBranchRateModel(p.getBranchRateModel());
	       setFrequencyModel(p.getFrequencyModel());
	}
}
