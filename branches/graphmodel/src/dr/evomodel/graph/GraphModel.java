package dr.evomodel.graph;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

import dr.evolution.alignment.SiteList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/*
 * A class to represent phylogenetic graphs where each node can have
 * up to two ancestors
 */
public class GraphModel extends TreeModel {

    public static final String GRAPH_MODEL = "graphModel";
	LinkedList<Node> freeNodes;	// a list of nodes for which storage exists 
	
    public GraphModel(Tree tree, PartitionModel partitionModel) {
        this(tree, false, partitionModel);
	}

	public GraphModel(String id, Tree tree, PartitionModel partitionModel) {
        this(tree, false, partitionModel);
        setId(id);
	}

    /* 
     * Creates a new GraphModel
     */
   public GraphModel(Tree tree, boolean copyAttributes, PartitionModel partitionModel) 
   {
	   // use the superconstructor but then convert all nodes
	   // from tree nodes to graph nodes
       super(tree);
       Node[] tmp = new Node[nodes.length];
       Node[] tmp2 = new Node[nodes.length];
       for(int i=0; i<tmp.length; i++)
       {
    	   tmp[i] = new Node();
    	   tmp2[i] = new Node();
    	   tmp[i].number = i;
    	   tmp2[i].number = i;
       }
       super.copyNodeStructure(tmp);
       for(int i=0; i<tmp.length; i++)
       {
    	   tmp[i].taxon = nodes[i].taxon;
       }
       nodes = storedNodes;
       super.copyNodeStructure(tmp2);
       for(int i=0; i<tmp.length; i++)
       {
    	   tmp2[i].taxon = nodes[i].taxon;
       }
       nodes = tmp;
       storedNodes = tmp2;
       root = nodes[root.number];
       
       // attach all siteRanges in the PartitionModel to each node
       for(int sr = 0; sr < partitionModel.getSiteRangeCount(); sr++){
           Partition range = partitionModel.getSiteRange(sr);
           for(int i=0; i<nodes.length; i++)
           {
        	   ((Node)nodes[i]).addObject(range);
           }
       }
   }
	
   /*
    * Add a new, unlinked node to the graph
    */
   public NodeRef newNode() {
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
       if(freeNodes.size()==0){
	       // need to expand storage to accommodate additional nodes
	       TreeModel.Node[] tmp = new TreeModel.Node[nodes.length*2];
	       TreeModel.Node[] tmp2 = new TreeModel.Node[storedNodes.length*2];
	       System.arraycopy(nodes, 0, tmp, 0, nodes.length);
	       System.arraycopy(storedNodes, 0, tmp2, 0, storedNodes.length);
	       for(int i=nodes.length; i<tmp.length; i++)
	       {
	    	   tmp[i] = new Node();
	    	   tmp[i].setNumber(i);
	           tmp[i].heightParameter = new Parameter.Default(0);
	           addVariable(tmp[i].heightParameter);
	    	   freeNodes.push((GraphModel.Node)tmp[i]);
	       }
	       for(int i=storedNodes.length; i<tmp2.length; i++)
	       {
	    	   tmp2[i] = new Node();
	    	   tmp2[i].setNumber(i);
	           tmp2[i].heightParameter = new Parameter.Default(0);
	           addVariable(tmp2[i].heightParameter);
	       }
	       nodes = tmp;
	       storedNodes = tmp2;
       }

       // simply return a node from the free list
       Node newNode = freeNodes.pop();
       internalNodeCount++;	// assume this is an internal node.  might not be true if there are partitions with subsets of taxa
       pushTreeChangedEvent(newNode);	// push a changed event onto the stack

       // add height, rate, and trait parameters
       // FIXME: do these parameters need to be created with default values?
       if(nhp!=null) nhp.addParameter(newNode.heightParameter);
       if(nrp!=null&&newNode.rateParameter!=null) nrp.addParameter(newNode.rateParameter);
       if(ntp!=null) {
           for (Map.Entry<String, Parameter> entry : newNode.getTraitMap().entrySet()) {
        	   ntp.addParameter(entry.getValue());
           }
       }
       
       return newNode;
   }

   /*
    * remove an unlinked node from the graph
    * @param node The node to remove
    */
   public void deleteNode(NodeRef node) {
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
       Node n = (Node) node;
       if(!n.hasNoChildren()||n.parent!=null||n.parent2!=null){
    	   throw new RuntimeException("Deleted node is linked to others!");
       }
       freeNodes.push(n);
       
       // remove from height, rate, and trait parameters
       if(nhp!=null) nhp.removeParameter(n.heightParameter);
       if(nrp!=null) nrp.removeParameter(n.rateParameter);
       if(ntp!=null) {
           for (Map.Entry<String, Parameter> entry : n.getTraitMap().entrySet()) {
        	   ntp.removeParameter(entry.getValue());
           }
       }
   }
   
   public void removePartition(NodeRef node, Partition range)
   {
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
	   Node n = (Node)node;
	   if(!n.hasObject(range)) throw new RuntimeException("Error, removing a nonexistant siterange!");
	   n.removeObject(range);
   }
   public void removePartitionFollowToRoot(NodeRef node, Partition range)
   {
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
	   // walk from node to root removing a site range
	   Node n = (Node)node;
	   while(n!=null){
		   if(n.hasObject(range)){
			   throw new RuntimeException("Error, removing a nonexistant partition!");
		   }
		   n.removeObject(range);
		   if(n.parent!=null && ((Node)n.parent).hasObject(range)){
			   n = (GraphModel.Node)n.parent;
		   }else if(n.parent2!=null && n.parent2.hasObject(range)){
			   n = (GraphModel.Node)n.parent2;
		   }else if(n.parent != null){
			   throw new RuntimeException("Error, no parent has relevant partition!");
		   }
	   }
   }

   public void addPartition(NodeRef node, Partition range)
   {
       if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
	   // walk from node to root removing a site range
	   Node n = (Node)node;
	   n.addObject(range);
	   pushTreeChangedEvent(n);
   }

   public void addSiteRangeFollowToRoot(NodeRef node, Partition range, Partition rangeToFollow)
   {
	   Node n = (Node)node;
	   while(n!=null){
		   n.addObject(range);
		   if(n.parent != null &&
				   ((Node)n.parent).hasObject(rangeToFollow)){
			   n = (Node)n.parent;
		   }else if(n.parent2 != null && 
				   n.parent2.hasObject(rangeToFollow)){
			   n = n.parent2;
		   }else{
			   throw new RuntimeException("Error, following a nonexistant siterange!");
		   }
	   }
   }
   /**
    * Uses TreeModel to copy the data, then links up the second parent
    */
   void copyNodeStructure(Node[] destination) {
	   super.copyNodeStructure(destination);

       for (int i = 0, n = nodes.length; i < n; i++) {
           Node node0 = (GraphModel.Node)nodes[i];
           Node node1 = destination[i];
           if (node0.parent2 != null) {
               node1.parent2 = (GraphModel.Node)storedNodes[node0.parent2.getNumber()];
           } else {
               node1.parent2 = null;
           }
       }
   }
   

   CompoundParameter nhp = null, nrp = null, ntp = null;
   public Parameter createNodeHeightsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {	   
	   CompoundParameter tmp = (CompoundParameter)super.createNodeHeightsParameter(rootNode, internalNodes, leafNodes);
	   if(internalNodes&&!rootNode) nhp = tmp;
	   return tmp;
   }
   public Parameter createNodeRatesParameter(double[] initialValues, boolean rootNode, boolean internalNodes, boolean leafNodes) {
	   nrp = (CompoundParameter)super.createNodeRatesParameter(initialValues, rootNode, internalNodes, leafNodes);
	   return nrp;
   }
   public Parameter createNodeTraitsParameter(String name, int dim, double[] initialValues,
           boolean rootNode, boolean internalNodes,
           boolean leafNodes, boolean firesTreeEvents) {
	   ntp = (CompoundParameter)createNodeTraitsParameter(name, dim, initialValues, rootNode, internalNodes, leafNodes, firesTreeEvents);
	   return ntp;
   }
   
   
   protected void handleModelChangedEvent(Model model, Object object, int index) 
   {
       // presumably a constituent partition has changed
   }
   
   /**
    * Converts all graph links to a string
    * @return a string
    */
   public String linkDump(){
	   StringBuilder sb = new StringBuilder();
	   for(int i=0; i<nodes.length; i++)
		   sb.append(((GraphModel.Node)nodes[i]).linksToString());
	   return sb.toString();
   }

    // **************************************************************
    // Private inner classes
    // **************************************************************

   	public class Node extends TreeModel.Node {

    	// a treeNode will have parent2 == null
    	// an argNode will have rightNode == null

    	public Node parent2 = null;	// an extra parent for recombinant nodes

    	protected HashSet<Object> objects;	// arbitrary objects tied to this node.  TODO: use the generic object mapper mentioned by Andrew
    	
    	
        public Node() {
        	super();
        	
        	objects = new HashSet<Object>();
        }

        public Node getChild(int n) {
        	return (Node)super.getChild(n);
        }

        /**
         * add new child node
         *
         * @param node new child node
         */
        public void addChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException("GraphModel.Nodes can only have 2 children");
            }
            if(node.parent==null){
            	node.parent = this;
            }else if(node.parent2==null){
            	node.parent2 = this;
            }else{
                throw new IllegalArgumentException("GraphModel.Nodes can only have 2 parents");
            }
        }

        /**
         * remove child
         *
         * @param node child to be removed
         */
        public Node removeChild(Node node) {
            if (leftChild == node) {
                leftChild = null;
            } else if (rightChild == node) {
                rightChild = null;
            } else {
                throw new IllegalArgumentException("Unknown child node");
            }
            if (node.parent == node) {
                node.parent = node.parent2;
                node.parent2 = null;
            } else if (node.parent2 == node) {
                node.parent2 = null;
            } else {
                throw new IllegalArgumentException("Unknown parent node");
            }
            return node;
        }

        /**
         * remove child
         *
         * @param n number of child to be removed
         */
        public Node removeChild(int n) {
            if (n == 0) {
                return removeChild((GraphModel.Node)leftChild);
            } else if (n == 1) {
                return removeChild((GraphModel.Node)rightChild);
            } else {
                throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
            }
        }
        
        public HashSet<Object> getObjects() {
        	return objects;
        }
        public boolean hasObject(Object o) {
        	return objects.contains(o);
        }
        public void addObject(Object o){
        	objects.add(o);
        }
        public void removeObject(Object o){
        	objects.remove(o);
        }


        public String toString() {
            return "node " + number + ", height=" + getHeight() + (taxon != null ? ": " + taxon.getId() : "");
        }
        
        public String linksToString() {
        	StringBuilder sb = new StringBuilder();
        	sb.append("node " + number);
        	if(parent!=null) sb.append(" parent1 " + parent.number);
        	if(parent2!=null) sb.append(" parent2 " + parent2.number);
        	if(leftChild!=null) sb.append(" leftChild " + leftChild.number);
        	if(rightChild!=null) sb.append(" rightChild " + rightChild.number);
        	return sb.toString();
        }
    }

}
