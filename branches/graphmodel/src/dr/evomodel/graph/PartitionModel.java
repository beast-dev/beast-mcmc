package dr.evomodel.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import dr.evolution.alignment.SiteList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/*
 * A PartitionModel stores one or more SiteLists (e.g. alignments) and a partitioning
 * over each SiteList into one or more contiguous segments.
 * Each Partition has an associated list of evolutionary models, such as SiteModels
 * or BranchRateModels or BranchSiteModels.
 * 
 */
public class PartitionModel extends AbstractModel {

    public static final String PARTITION_MODEL = "partitionModel";
	
    protected List<SiteList> siteLists;	// the set of one or more siteLists on which this partition exists

	protected Partition[] partitions;
	protected Partition[] storedPartitions;
	protected LinkedList<Partition> freePartitions;

	protected HashMap<Partition,List<Model>> modelsOnPartition;
	
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
		partitions = new Partition[siteLists.size()];
		storedPartitions = new Partition[siteLists.size()];
		for(int i=0; i<siteLists.size(); i++)
		{
			partitions[i] = new Partition(siteLists.get(i));
			partitions[i].setNumber(i);
			storedPartitions[i] = new Partition(siteLists.get(i));
			storedPartitions[i].setNumber(i);
		}
		modelsOnPartition = new HashMap<Partition, List<Model>>();
	}
	
    public void pushPartitionChangedEvent(Partition siteRange, int left, int right) {
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
        // no partitions overlap

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
	Partition newPartition(Partition siteRange)
	{
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
       if(freePartitions.size()==0){
	       // need to expand storage to accommodate additional nodes
	       Partition[] tmp = new Partition[partitions.length*2];
	       Partition[] tmp2 = new Partition[storedPartitions.length*2];
	       System.arraycopy(partitions, 0, tmp, 0, partitions.length);
	       System.arraycopy(storedPartitions, 0, tmp2, 0, storedPartitions.length);
	       for(int i=partitions.length; i<tmp.length; i++)
	       {
	    	   tmp[i] = new Partition();
	    	   tmp[i].setNumber(i);
	    	   freePartitions.push(tmp[i]);
	       }
	       for(int i=storedPartitions.length; i<tmp2.length; i++)
	       {
	    	   tmp2[i] = new Partition();
	    	   tmp2[i].setNumber(i);
	       }
	       partitions = tmp;
	       storedPartitions = tmp2;
       }

       // get a new SiteRange and copy the values of the provided one
       Partition newSR = freePartitions.pop();
       newSR.setLeftSite(siteRange.getLeftSite());
       newSR.setRightSite(siteRange.getRightSite());
       newSR.setSiteList(siteRange.getSiteList());

       // add a model list for it
       ArrayList<Model> al = new ArrayList<Model>();
       modelsOnPartition.put(newSR, al);

       return newSR;
	}
	
	void changeRange(Partition siteRange, int newLeft, int newRight)
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
	
	public void removeSiteRange(Partition siteRange){
		freePartitions.push(siteRange);
		modelsOnPartition.remove(siteRange);
	}
	
	public Partition getSiteRange(int i){
		return partitions[i];
	}
	public int getSiteRangeCount(){
		return partitions.length;
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
	
	public List<Model> getModelsOnPartition(Partition siteRange)
	{
		return modelsOnPartition.get(siteRange);
	}
	
	void addModelToPartition(Partition siteRange, Model model){
		List<Model> l = modelsOnPartition.get(siteRange);
		l.add(model);
	}
	
	void removeModelFromPartition(Partition siteRange, Model model){
		List<Model> l = modelsOnPartition.get(siteRange);
		l.remove(model);
	}

	@Override
	protected void restoreState() {
		partitions = storedPartitions;
	}

	@Override
	protected void storeState() {
		for(int i = 0; i<partitions.length; i++){
			storedPartitions[i].setLeftSite(partitions[i].getLeftSite());
			storedPartitions[i].setRightSite(partitions[i].getRightSite());
			storedPartitions[i].setSiteList(partitions[i].getSiteList());
		}
	}

    public class PartitionChangedEvent {
        final Partition partition;
        final int newSectionLeft;
        final int newSectionRight;

        public PartitionChangedEvent() {
            this(null, -1, -1);
        }

        public PartitionChangedEvent(Partition partition) {
            this(partition, -1, -1);
        }

        public PartitionChangedEvent(Partition partition, int newSectionLeft, int newSectionRight) {
            this.partition = partition;
            this.newSectionLeft = newSectionLeft;
            this.newSectionRight = newSectionRight;
        }

        public Partition getPartition() {
            return partition;
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
