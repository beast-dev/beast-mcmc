/*
 * ARGTree.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.arg;

import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.arg.ARGModel.Node;
import dr.util.Attributable;

import java.util.*;

public class ARGTree implements Tree {

	//NodeRef root;

	protected Taxon[] taxaList;
	protected int taxaCount;

	private final Node initialRoot;

    public ARGModel argModel;
//    private Map<NodeRef,Integer> mapARGNodesToInts;

    private final Map<NodeRef,NodeRef> mapARGNodesToTreeNodes;

//    public Map<NodeRef,Integer> getMappingInts() { return mapARGNodesToInts; }

//    public Map<ARGModel.Node, NodeRef> getMappingNodes() { return mapARGNodesToTreeNodes; }

	private int partition = -9;

//	/**
//	 * Constructor to represent complete ARG as a tree
//	 *
//	 * @param arg
//	 */
//	public ARGTree(ARGModel arg) {
//              this.argModel = arg;
//              mapNodesARGToTree = new HashMap<Node,Node>(arg.getNodeCount());
//              root = arg.new Node((Node) arg.getRoot());
//	}


//	/**
//	 * Constructor for specific partition tree
//	 *
//	 * @param arg
//	 * @param partition
//	 */

//	ArrayList<Node> nodeList;

	public boolean wasRootTrimmed() {
		return (root != initialRoot);
	}


	public String toGraphString() {
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

	public ARGTree(ARGModel arg, int partition) {
              this.argModel = arg;
              mapARGNodesToTreeNodes = new HashMap<NodeRef,NodeRef>(arg.getNodeCount());

		this.partition = partition;
		ARGModel.Node node = arg.new Node(((Node) arg.getRoot()), partition);
		initialRoot = node;

		int j = arg.externalNodeCount;
		node.stripOutDeadEnds();
		root = node.stripOutSingleChildNodes(node);
		node = root;
		nodeCount = 2 * j - 1;
		externalNodeCount = j;
		internalNodeCount = j - 1;
		nodes = new Node[nodeCount];

		do {
			node = (Node) TreeUtils.postorderSuccessor(this, node);
			if (node.isExternal()) {
                                  // keep same order as ARG, so do not need to reload tipStates/Partials
                                  nodes[node.number] = node;
                                  mapARGNodesToTreeNodes.put(node.mirrorNode,node);
			} else {
                                  // Reorder in new post-order succession

				nodes[j] = node;
                                        node.number = j;
				j++;
                                  mapARGNodesToTreeNodes.put(node.mirrorNode, node);
			}
		} while (node != root);
        hasRates = false;
    }


    public Map<NodeRef,NodeRef> getMapping() {
        return mapARGNodesToTreeNodes;
    }

//    public Map<NodeRef,Integer> getMapARGNodesToInts() {
//        // Only need to map internal nodes
//        Map<NodeRef,Integer> map = new HashMap<NodeRef,Integer>(getInternalNodeCount());
//        for(int i=0; i<getInternalNodeCount(); i++) {
//            ARGModel.Node node = (ARGModel.Node) getInternalNode(i);
//            map.put(node.mirrorNode,node.number);
//        }
//        return map;
//    }
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

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index ++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
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

	private Attributable.AttributeHelper treeAttributes = null;

	/**
	 * Sets an named attribute for this object.
	 *
	 * @param name  the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setAttribute(String name, Object value) {
		if (treeAttributes == null)
			treeAttributes = new Attributable.AttributeHelper();
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
	public Iterator<String> getAttributeNames() {
		if (treeAttributes == null)
			return null;
		else
			return treeAttributes.getAttributeNames();
	}

	/**
	 * @return a string containing a newick representation of the tree
	 */
	public final String getNewick() {
		return TreeUtils.newick(this);
	}

	public final String getUniqueNewick(){
		return TreeUtils.uniqueNewick(this,this.getRoot());
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

	private final boolean hasRates;
	private final boolean hasTraits = false;


}
