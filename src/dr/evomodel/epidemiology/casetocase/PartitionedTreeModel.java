package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.*;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import java.util.Set;

/**
 * TreeModel plus partition information
 *
 * todo a lot of methods should eventually move here
 */
public class PartitionedTreeModel extends TreeModel {

    private BranchMapModel branchMap;
    public final static String PARTITIONED_TREE_MODEL = "partitionedTreeModel";
    Set<NodeRef> partitionsQueue = new HashSet<NodeRef>();

    public PartitionedTreeModel(String id, Tree tree, BranchMapModel branchMap){
        super(id, tree);
        this.branchMap = branchMap;
        addModel(branchMap);
    }

    public PartitionedTreeModel(String id, Tree tree){
        super(id, tree);
        branchMap = new BranchMapModel(this);
        addModel(branchMap);
    }

    public PartitionedTreeModel(TreeModel treeModel, BranchMapModel branchMap){
        this(PARTITIONED_TREE_MODEL, treeModel, branchMap);
    }

    public PartitionedTreeModel(TreeModel treeModel){
        this(PARTITIONED_TREE_MODEL, treeModel);
    }

    public void partitionsChangingAlert(HashSet<AbstractCase> casesToRecalculate){
        listenerHelper.fireModelChanged(this, new PartitionsChangedEvent(casesToRecalculate));
    }

    public void partitionChangingAlert(AbstractCase caseToRecalculate){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        out.add(caseToRecalculate);

        partitionsChangingAlert(out);
    }

    public void universalAlert(){
        HashSet<AbstractCase> allCases = new HashSet<AbstractCase>(Arrays.asList(branchMap.getArrayCopy()));
        partitionsChangingAlert(allCases);
    }

    public BranchMapModel getBranchMap(){
        return branchMap;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        // should only happen when the branch map changes

        if(object instanceof ArrayList){
            HashSet<AbstractCase> changedPartitions = new HashSet<AbstractCase>();

            for(int i=0; i<((ArrayList) object).size(); i++){
                BranchMapModel.BranchMapChangedEvent event
                        =  (BranchMapModel.BranchMapChangedEvent)((ArrayList) object).get(i);

                changedPartitions.add(event.getOldCase());
                changedPartitions.add(event.getNewCase());

                NodeRef node = getNode(event.getNodeToRecalculate());
                NodeRef parent = getParent(node);

                changedPartitions.add(branchMap.get(parent.getNumber()));


            }
            partitionsChangingAlert(changedPartitions);


        } else {
            throw new RuntimeException("Unanticipated model changed event from BranchMapModel");
        }

    }

    public class PartitionsChangedEvent{

        private final HashSet<AbstractCase> casesToRecalculate;

        public PartitionsChangedEvent(HashSet<AbstractCase> casesToRecalculate){
            this.casesToRecalculate = casesToRecalculate;
        }

        public HashSet<AbstractCase> getCasesToRecalculate(){
            return casesToRecalculate;
        }
    }


    public void pushNodePartitionsChangedEvent(NodeRef node){
        int nodeNumber = node.getNumber();

        if(!inTreeEdit()){
            partitionChangingAlert(branchMap.get(nodeNumber));
        } else {
            partitionsQueue.add(node);
        }
    }


    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
            //this signals events where values in all dimensions of a parameter is changed.
            universalAlert();
        } else {

            final NodeRef node = getNodeOfParameter((Parameter) variable);

            partitionsChangingAlert(adjacentPartitions(node));
        }
    }

    public HashSet<AbstractCase> adjacentPartitions(NodeRef node){
        HashSet<AbstractCase> changedCases = new HashSet<AbstractCase>();
        ArrayList<NodeRef> affectedNodes = new ArrayList<NodeRef>();

        affectedNodes.add(node);
        affectedNodes.add(getParent(node));
        affectedNodes.add(getChild(node, 0));
        affectedNodes.add(getChild(node, 1));

        for(NodeRef aNode : affectedNodes){
            changedCases.add(branchMap.get(aNode.getNumber()));
        }

        return changedCases;
    }

    private void flushQueue(){
        if(!inTreeEdit()){
            throw new RuntimeException("Wait until you've finished editing the tree before flushing the partition" +
                    "queue");
        }

        for(NodeRef node : partitionsQueue){
            AbstractCase nodePartition = branchMap.get(node.getNumber());
            partitionChangingAlert(nodePartition);
            NodeRef parent = getParent(node);
            if(parent!=null && branchMap.get(node.getNumber())!=branchMap.get(parent.getNumber())){
                partitionChangingAlert(branchMap.get(parent.getNumber()));
            }
        }

        partitionsQueue.clear();

    }

    public void addChild(NodeRef p, NodeRef c) {

        partitionChangingAlert(branchMap.get(p.getNumber()));
        partitionChangingAlert(branchMap.get(c.getNumber()));

        super.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        partitionChangingAlert(branchMap.get(p.getNumber()));
        partitionChangingAlert(branchMap.get(c.getNumber()));

        super.removeChild(p, c);

    }

    public void setNodeHeight(NodeRef n, double height) {

        partitionsChangingAlert(adjacentPartitions(n));

        super.setNodeHeight(n, height);
    }



    // anything you do to the partitions must be finished before you call this.

    public void endTreeEdit() {
        super.endTreeEdit();

        // todo in the end, want to check the tree partitions are sane here before flushing the queue

        flushQueue();
    }

}
