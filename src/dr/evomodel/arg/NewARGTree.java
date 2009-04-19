package dr.evomodel.arg;

import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.arg.ARGModel.Node;

import java.util.ArrayList;
import java.util.Iterator;

public class NewARGTree implements Tree {

    //NodeRef root;

    protected Taxon[] taxaList;
    protected int taxaCount;

    private Node initialRoot;

    private int partition = -9;

    /**
     * Constructor to represent complete ARG as a tree
     *
     * @param arg
     */
    public NewARGTree(ARGModel arg) {
        root = arg.new Node((Node) arg.getRoot());
//		  root = arg.new Node( (Node)arg.getRoot() );
        // int cnt = countReassortmentNodes(root);
    }


    /**
     * Constructor for specific partition tree
     *
     * @param arg
     * @param partition
     */

    ArrayList<Node> nodeList;

    public boolean wasRootTrimmed() {
        return (root != initialRoot);
//		return false;
    }


    public String toGraphString() {
        //if( true )
        //	  return null;
//		  int number = 1;
//		  for( Node node : nodeList )
//			  node.number = number++;
        StringBuffer sb = new StringBuffer();
        for (Node node : nodes) {
            sb.append(node.number);
            if (node.leftParent != null)
                sb.append(" " + node.leftParent.number);
            else
                sb.append(" 0");
            if (node.rightParent != null)
                sb.append(" " + node.rightParent.number);
            else
                sb.append(" 0");
            if (node.leftChild != null)
                sb.append(" " + node.leftChild.number);
            else
                sb.append(" 0");
            if (node.rightChild != null)
                sb.append(" " + node.rightChild.number);
            else
                sb.append(" 0");
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
            sb.append("\n");
        }
        sb.append("Root = " + ((Node) getRoot()).number + "\n");
        return new String(sb);
    }

    public NewARGTree(ARGModel arg, int partition) {

        this.partition = partition;


        nodeList = new ArrayList<Node>();
        Node node = arg.new Node(((Node) arg.getRoot()), partition, nodeList);
        initialRoot = node;


        int i = 0;
        int j = arg.externalNodeCount;

        //root = node;
        //System.err.println("Null rights 1 = "+checkForNullRights(node));
//		System.err.println("ARG->TREE1\n"+this.toGraphString());
        node.stripOutDeadEnds();
//		System.err.println("ARG->TREE2\n"+this.toGraphString());
        root = node.stripOutSingleChildNodes(node);
//		System.err.println("root == node? "+(root == node));
        node = root;
        //System.err.println("Null rights 2 - "+checkForNullRights(node));
//		System.err.println("ARG->TREE3\n"+this.toGraphString());
        //root = node;
        //if( save != root ) {
        //	System.err.println("Rerooted: "+Tree.Utils.uniqueNewick(this, node));
        //	System.exit(-1);
        //}

        //Tree.Utils.uniqueNewick(this, root));
        //	System.exit(-1);
//		System.err.println("root.number = "+root.number);
        nodeCount = 2 * j - 1;
        externalNodeCount = j;
        internalNodeCount = j - 1;
        nodes = new Node[nodeCount];


        do {
            node = (Node) Utils.postorderSuccessor(this, node);
            //System.err.println("Ordering: "+Tree.Utils.uniqueNewick(this, node));
            if (node.isExternal()) {
                node.number = i;

                nodes[i] = node;
                // storedNodes[i] = new Node();
                // storedNodes[i].taxon = node.taxon;
                // storedNodes[i].number = i;

                i++;
            } else {
                node.number = j;

                nodes[j] = node;
                // storedNodes[j] = new Node();
                // storedNodes[j].number = j;

                j++;
            }
        } while (node != root);

//		System.err.println("After: "+this.toString());
        //System.err.println("Finished post-order on ARGTree");
    }

    public boolean checkForNullRights(Node node) {
        return node.checkForNullRights();
    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    /**
     * Return the units that this tree is expressed in.
     */

    public final Type getUnits() {
        return units;
    }

    public void setUnits(Type units) {
        this.units = units;

    }

    /**
     * Sets the units that this tree is expressed in.
     */


    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree.
     */
    public final int getNodeCount() {
        return nodeCount;
    }

    public final boolean hasNodeHeights() {
        return true;
    }

    public final double getNodeHeight(NodeRef node) {

        //System.err.println(Tree.Utils.uniqueNewick(this, node));
        //((Node)node))

        return ((Node) node).getHeight();
    }

    public final double getNodeHeightUpper(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getUpperLimit(0);
    }

    public final double getNodeHeightLower(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getLowerLimit(0);
    }


    /**
     * @param node
     * @return the rate parameter associated with this node.
     */
    public final double getNodeRate(NodeRef node) {


        if (!hasRates) {
            return 1.0;
        }

        return ((Node) node).getRate(partition);
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
    }

    public double getNodeTrait(NodeRef node) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getTrait();
    }

    public final Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public final boolean isExternal(NodeRef node) {
        return ((Node) node).isExternal();
    }

    public final boolean isRoot(NodeRef node) {
        return (node == root);
    }

    public final int getChildCount(NodeRef node) {
        //System.err.println("Cn for "+((Node)node).number);
        return ((Node) node).getChildCount();
    }

    public final NodeRef getChild(NodeRef node, int i) {
        return ((Node) node).getChild(i);
    }

    public final NodeRef getParent(NodeRef node) {
        //System.err.println("Gimme!");
        return ((Node) node).leftParent;
    }

    public final boolean hasBranchLengths() {
        return true;
    }

    public final double getBranchLength(NodeRef node) {
        NodeRef parent = getParent(node);
        if (parent == null) {
            return 0.0;
        }

        return getNodeHeight(parent) - getNodeHeight(node);
    }

    public final NodeRef getExternalNode(int i) {
        return nodes[i];
    }

    public final NodeRef getInternalNode(int i) {
        return nodes[i + externalNodeCount];
    }

    public final NodeRef getNode(int i) {
        return nodes[i];
    }

    /**
     * Returns the number of external nodes.
     */
    public final int getExternalNodeCount() {
        return externalNodeCount;
    }

    /**
     * Returns the ith internal node.
     */
    public final int getInternalNodeCount() {
        return internalNodeCount;
    }

    /**
     * Returns the root node of this tree.
     */
    public final NodeRef getRoot() {
        return root;
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        return getExternalNodeCount();
    }

    /**
     * @return the ith taxon in the list.
     */
//    public Taxon getTaxon(int taxonIndex) {
//        return ((Node)getExternalNode(taxonIndex)).taxon;
//    }

//    public Taxon getTaxon(int taxonIndex) {
//    	if( taxonIndex >= taxaCount )
//    		return null;
//        return taxaList[taxonIndex];
//    }
    public Taxon getTaxon(int taxonIndex) {
        return ((Node) getExternalNode(taxonIndex)).taxon;
    }


    /**
     * @return the ID of the taxon of the ith external node. If it doesn't have
     *         a taxon, returns the ID of the node itself.
     */
    public String getTaxonId(int taxonIndex) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getId();
        } else {
            return null;
        }
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon) return i;
        }
        return -1;
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the taxon of the given
     *         external node. If the node doesn't have a taxon then the nodes own attribute
     *         is returned.
     */
    public final Object getTaxonAttribute(int taxonIndex, String name) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getAttribute(name);
        }
        return null;
    }

    // **************************************************************
    // MutableTaxonList IMPLEMENTATION
    // **************************************************************

    public int addTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public void setTaxonId(int taxonIndex, String id) {
        throw new IllegalArgumentException("Cannot set taxon id in a TreeModel");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new IllegalArgumentException("Cannot set taxon attribute in a TreeModel");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
    } // Do nothing at the moment

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
    } // Do nothing at the moment

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    protected String id = null;

    /**
     * @return the id.
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the id.
     */
    public void setId(String id) {
        this.id = id;
    }

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    private AttributeHelper treeAttributes = null;

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (treeAttributes == null)
            treeAttributes = new AttributeHelper();
        treeAttributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator getAttributeNames() {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttributeNames();
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public final String getNewick() {
        return Utils.newick(this);
    }

    public final String getUniqueNewick() {
        return Utils.uniqueNewick(this, this.getRoot());
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public String toString() {
        return getNewick();
    }

    public Tree getCopy() {
        throw new UnsupportedOperationException("please don't call this function");
    }

    // ***********************************************************************
    // Private members
    // ***********************************************************************


    /**
     * root node
     */
    protected Node root = null;
    protected int storedRootNumber;

    /**
     * list of internal nodes (including root)
     */
    protected Node[] nodes = null;
    protected Node[] storedNodes = null;

    /**
     * number of nodes (including root and tips)
     */
    protected int nodeCount;

    /**
     * number of external nodes
     */
    protected int externalNodeCount;

    /**
     * number of internal nodes (including root)
     */
    protected int internalNodeCount;


    /**
     * holds the units of the trees branches.
     */
//	private int units = SUBSTITUTIONS;
    private Type units;

    protected boolean inEdit = false;

    private boolean hasRates = false;
    private boolean hasTraits = false;


}