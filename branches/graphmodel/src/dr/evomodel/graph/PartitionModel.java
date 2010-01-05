package dr.evomodel.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import dr.evolution.alignment.SiteList;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeModel.Node;
import dr.evomodel.tree.TreeModel.TreeChangedEvent;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/*
 * how to do this?
 * should PartitionModel keep a set of SiteRange
 * and a separate list of SiteModel?
 * Then each SiteRange would have a reference to
 * the corresponding SiteModel and SiteList?
 * would have modifiers to add and remove SiteRanges
 * change the SiteModel of a SiteRange, or change the
 * range of a SiteRange
 * When a range of a SiteRange changes, one pushes a SiteRangeChanged event
 * When a new SiteRange emerges, one pushes a NewSiteRange event
 * When a constituent SiteModel changes, the event is propagated
 * When a SiteRange's SiteModel is changed, a SiteRangeSiteModelChanged event
 */
public class PartitionModel extends AbstractModel {

    public static final String PARTITION_MODEL = "partitionModel";
	
    protected List<SiteList> siteLists;	// the set of one or more siteLists on which this partition exists

	protected SiteRange[] siteRanges;
	protected SiteRange[] storedSiteRanges;
	protected LinkedList<SiteRange> freeSiteRanges;

	protected HashMap<SiteRange,List<Model>> siteRangeModels;
	
    protected boolean inEdit = false;
    protected final List<PartitionChangedEvent> partitionChangedEvents = new ArrayList<PartitionChangedEvent>();

    /**
     * Construct a Partition Model over one or more SiteList objects (alignments)
     * @param siteLists
     */
	public PartitionModel(List<SiteList> siteLists)
	{
		super(PARTITION_MODEL);
		
		this.siteLists = new ArrayList<SiteList>(siteLists.size());
		
		this.siteLists.addAll(siteLists);
		siteRanges = new SiteRange[siteLists.size()];
		storedSiteRanges = new SiteRange[siteLists.size()];
		for(int i=0; i<siteLists.size(); i++)
		{
			siteRanges[i] = new SiteRange(siteLists.get(i));
			siteRanges[i].setNumber(i);
			storedSiteRanges[i] = new SiteRange(siteLists.get(i));
			storedSiteRanges[i].setNumber(i);
		}
		siteRangeModels = new HashMap<SiteRange, List<Model>>();
	}
	
    public void pushPartitionChangedEvent(SiteRange siteRange, int left, int right) {
    	PartitionChangedEvent pce = new PartitionChangedEvent(siteRange, left, right);
    	pushPartitionChangedEvent(pce);
    }
    public void pushPartitionChangedEvent(PartitionChangedEvent event) {
        if (inEdit) {
            partitionChangedEvents.add(event);
        } else {
            listenerHelper.fireModelChanged(this, event);
        }
    }

    public boolean beginEdit() {
        if (inEdit) throw new RuntimeException("Alreading in edit transaction mode!");
        inEdit = true;
        return false;
    }

    public void endEdit() {
        if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");
        inEdit = false;

        // TODO:
        // check that the partition over sites is complete and
        // no siteranges overlap

        for (PartitionChangedEvent partitionChangedEvent : partitionChangedEvents) {
            listenerHelper.fireModelChanged(this, partitionChangedEvent);
        }
        partitionChangedEvents.clear();
    }
	
	/*
	 * Allocate a new SiteRange by copying an existing SiteRange
	 * The new SiteRange will overlap, so it is necessary to 
	 * change the boundaries of new and/or old SiteRanges
	 */
	SiteRange newSiteRange(SiteRange siteRange)
	{
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
       if(freeSiteRanges.size()==0){
	       // need to expand storage to accommodate additional nodes
	       SiteRange[] tmp = new SiteRange[siteRanges.length*2];
	       SiteRange[] tmp2 = new SiteRange[storedSiteRanges.length*2];
	       System.arraycopy(siteRanges, 0, tmp, 0, siteRanges.length);
	       System.arraycopy(storedSiteRanges, 0, tmp2, 0, storedSiteRanges.length);
	       for(int i=siteRanges.length; i<tmp.length; i++)
	       {
	    	   tmp[i] = new SiteRange();
	    	   tmp[i].setNumber(i);
	    	   freeSiteRanges.push(tmp[i]);
	       }
	       for(int i=storedSiteRanges.length; i<tmp2.length; i++)
	       {
	    	   tmp2[i] = new SiteRange();
	    	   tmp2[i].setNumber(i);
	       }
	       siteRanges = tmp;
	       storedSiteRanges = tmp2;
       }

       // get a new SiteRange and copy the values of the provided one
       SiteRange newSR = freeSiteRanges.pop();
       newSR.setLeftSite(siteRange.getLeftSite());
       newSR.setRightSite(siteRange.getRightSite());
       newSR.setSiteList(siteRange.getSiteList());

       // add a model list for it
       ArrayList<Model> al = new ArrayList<Model>();
       siteRangeModels.put(newSR, al);

       return newSR;
	}
	
	void changeRange(SiteRange siteRange, int newLeft, int newRight)
	{
	       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
		// check whether the range has expanded or shrunk
		// if expanded, we need to compute new likelihoods
		if(newLeft < siteRange.getLeftSite()){
			int minnow = newRight < siteRange.getLeftSite() ? newRight : siteRange.getLeftSite();
			pushPartitionChangedEvent(siteRange, newLeft, minnow);
		}
		if(newRight > siteRange.getRightSite()){
			int maxxow = newLeft > siteRange.getRightSite() ? newLeft : siteRange.getRightSite();
			pushPartitionChangedEvent(siteRange, maxxow, siteRange.getRightSite());
		}
		siteRange.setLeftSite(newLeft);
		siteRange.setRightSite(newRight);
	}
	
	public void removeSiteRange(SiteRange siteRange){
		freeSiteRanges.push(siteRange);
		siteRangeModels.remove(siteRange);
	}
	
	public SiteRange getSiteRange(int i){
		return siteRanges[i];
	}
	public int getSiteRangeCount(){
		return siteRanges.length;
	}
	
	public SiteList getSiteList(int i){
		return siteLists.get(i);
	}
	public int getSiteListCount(){
		return siteLists.size();
	}

	@Override
	protected void acceptState() {
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub

	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub

	}
	
	public List<Model> getSiteRangeModels(SiteRange siteRange)
	{
		return siteRangeModels.get(siteRange);
	}
	
	void addSiteRangeModel(SiteRange siteRange, SiteModel model){
		List<Model> l = siteRangeModels.get(siteRange);
		l.add(model);
	}
	
	void removeSiteRangeModel(SiteRange siteRange, SiteModel model){
		List<Model> l = siteRangeModels.get(siteRange);
		l.remove(model);
	}

	@Override
	protected void restoreState() {
		siteRanges = storedSiteRanges;
	}

	@Override
	protected void storeState() {
		for(int i = 0; i<siteRanges.length; i++){
			storedSiteRanges[i].setLeftSite(siteRanges[i].getLeftSite());
			storedSiteRanges[i].setRightSite(siteRanges[i].getRightSite());
			storedSiteRanges[i].setSiteList(siteRanges[i].getSiteList());
		}
	}

    public class PartitionChangedEvent {
        final SiteRange siteRange;
        final int newSectionLeft;
        final int newSectionRight;

        public PartitionChangedEvent() {
            this(null, -1, -1);
        }

        public PartitionChangedEvent(SiteRange siteRange) {
            this(siteRange, -1, -1);
        }

        public PartitionChangedEvent(SiteRange siteRange, int newSectionLeft, int newSectionRight) {
            this.siteRange = siteRange;
            this.newSectionLeft = newSectionLeft;
            this.newSectionRight = newSectionRight;
        }

        public SiteRange getSiteRange() {
            return siteRange;
        }
        public boolean hasNewSection() {
            return newSectionRight != -1;
        }
        public int getNewSectionLeft() {
        	return newSectionLeft;
        }
        public int getNewSectionRight() {
        	return newSectionLeft;
        }
    }
}
