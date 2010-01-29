package dr.evomodel.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import dr.evolution.alignment.SiteList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
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
	protected HashMap<Partition,List<Model>> storedModelsOnPartition;
	
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
		modelsOnPartition = new HashMap<Partition, List<Model>>();
		for(int i=0; i<siteLists.size(); i++)
		{
			partitions[i] = new Partition(siteLists.get(i));
			partitions[i].setNumber(i);
			storedPartitions[i] = new Partition(siteLists.get(i));
			storedPartitions[i].setNumber(i);
			modelsOnPartition.put(partitions[i], new ArrayList<Model>());
			modelsOnPartition.put(storedPartitions[i], new ArrayList<Model>());
		}
	}
	
    public void pushPartitionChangedEvent(Partition partition, int left, int right) {
    	PartitionChangedEvent pce = new PartitionChangedEvent(partition, left, right);
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
	 * Allocate a new Partition by copying an existing Partition
	 * The new Partition will overlap, so it is necessary to 
	 * change the boundaries of new and/or old Partitions
	 */
	Partition newPartition(Partition partition)
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

       // get a new Partition and copy the values of the provided one
       Partition newSR = freePartitions.pop();
       newSR.copyPartition(partition);

       // add a model list for it
       ArrayList<Model> al = new ArrayList<Model>();
       modelsOnPartition.put(newSR, al);

       return newSR;
	}
	
	void changeRange(Partition partition, int newLeft, int newRight)
	{
	       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
		// check whether the range has expanded or shrunk
		// if expanded, we need to compute new likelihoods
		if(newLeft < partition.getLeftSite()){
			int minnow = newRight < partition.getLeftSite() ? newRight : partition.getLeftSite();
			pushPartitionChangedEvent(partition, newLeft, minnow);
		}
		if(newRight > partition.getRightSite()){
			int maxxow = newLeft > partition.getRightSite() ? newLeft : partition.getRightSite();
			pushPartitionChangedEvent(partition, maxxow, partition.getRightSite());
		}
		partition.setLeftSite(newLeft);
		partition.setRightSite(newRight);
	}
	
	public void removePartition(Partition partition){
		freePartitions.push(partition);
		modelsOnPartition.remove(partition);
	}
	
	public Partition getPartition(int i){
		return partitions[i];
	}
	public int getPartitionCount(){
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
	
	public List<Model> getModelsOnPartition(Partition partition)
	{
		return modelsOnPartition.get(partition);
	}
	
	public void addModelToPartition(Partition partition, Model model){
		List<Model> l = modelsOnPartition.get(partition);
		if(model instanceof SiteModel)
			partition.setSiteModel((SiteModel)model);
		if(model instanceof BranchRateModel)
			partition.setBranchRateModel((BranchRateModel)model);
		l.add(model);
	}
	
	public void removeModelFromPartition(Partition partition, Model model){
		List<Model> l = modelsOnPartition.get(partition);
		l.remove(model);
	}

	@Override
	protected void restoreState() {
		partitions = storedPartitions;
		modelsOnPartition = storedModelsOnPartition;
	}

	@Override
	protected void storeState() {
		for(int i = 0; i<partitions.length; i++){
			storedPartitions[i].copyPartition(partitions[i]);
		}
		storedModelsOnPartition = modelsOnPartition;
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
