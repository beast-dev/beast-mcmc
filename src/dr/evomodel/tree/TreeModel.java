/*
 * TreeModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.*;
import dr.util.Attributable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.*;

/**
 * A model component for trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeModel.java,v 1.129 2006/01/05 17:55:47 rambaut Exp $
 */
public class TreeModel extends AbstractModel implements MutableTree {

	//
	// Public stuff
	//

	public static final String TREE_MODEL = "treeModel";

	public TreeModel(Tree tree) {

		super(TREE_MODEL);

		// get a rooted version of the tree to clone
		FlexibleTree binaryTree = new FlexibleTree(tree);
		binaryTree.resolveTree();

		// clone the node structure (this will create the individual parameters
		Node node = new Node(binaryTree, binaryTree.getRoot());

		internalNodeCount = binaryTree.getInternalNodeCount();
		externalNodeCount = binaryTree.getExternalNodeCount();

		nodeCount = internalNodeCount + externalNodeCount;

		nodes = new Node[nodeCount];
		storedNodes = new Node[nodeCount];

		int i = 0;
		int j = externalNodeCount;

		root = node;

		do {
			node = (Node) Tree.Utils.postorderSuccessor(this, node);

			if (node.isExternal()) {
				node.number = i;

				nodes[i] = node;
				storedNodes[i] = new Node();
				storedNodes[i].taxon = node.taxon;
				storedNodes[i].number = i;

				i++;
			} else {
				node.number = j;

				nodes[j] = node;
				storedNodes[j] = new Node();
				storedNodes[j].number = j;

				j++;
			}
		} while (node != root);

	}

	/* New constructor that copies the attributes of Tree tree into the new TreeModel
	 * Useful for constructing a TreeModel from a NEXUS file entry
	 */

	public TreeModel(Tree tree, boolean copyAttributes) {

		super(TREE_MODEL);

		// get a rooted version of the tree to clone
		FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
		binaryTree.resolveTree();

		// clone the node structure (this will create the individual parameters
		Node node = new Node(binaryTree, binaryTree.getRoot());

		internalNodeCount = binaryTree.getInternalNodeCount();
		externalNodeCount = binaryTree.getExternalNodeCount();

		nodeCount = internalNodeCount + externalNodeCount;

		nodes = new Node[nodeCount];
		storedNodes = new Node[nodeCount];

		int i = 0;
		int j = externalNodeCount;

		root = node;

		do {
			node = (Node) Tree.Utils.postorderSuccessor(this, node);

			if (node.isExternal()) {
				node.number = i;

				nodes[i] = node;
				storedNodes[i] = new Node();
				storedNodes[i].taxon = node.taxon;
				storedNodes[i].number = i;

				i++;
			} else {
				node.number = j;

				nodes[j] = node;
				storedNodes[j] = new Node();
				storedNodes[j].number = j;

				j++;
			}
		} while (node != root);

	}


	void setupHeightBounds() {

		for (int i = 0; i < nodeCount; i++) {
			nodes[i].setupHeightBounds();
		}
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent() {
		pushTreeChangedEvent(new TreeChangedEvent());
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(NodeRef nodeRef) {
		pushTreeChangedEvent(new TreeChangedEvent((Node) nodeRef));
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(Node node, Parameter parameter, int index) {
		pushTreeChangedEvent(new TreeChangedEvent(node, parameter, index));
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(TreeChangedEvent event) {
		if (inEdit) {
			treeChangedEvents.add(event);
		} else {
			listenerHelper.fireModelChanged(this, event);
		}
	}


	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no submodels so nothing to do
	}

	/**
	 * Called when a parameter changes.
	 */
	public void handleParameterChangedEvent(Parameter parameter, int index) {
		Node node = getNodeOfParameter(parameter);
		pushTreeChangedEvent(node, parameter, index);
	}

	private List<TreeChangedEvent> treeChangedEvents = new ArrayList<TreeChangedEvent>();

	public boolean hasRates() {
		return hasRates;
	}

	public class TreeChangedEvent {

		Node node;
		Parameter parameter;
		int index;

		public TreeChangedEvent() {
			this(null, null, -1);
		}

		public TreeChangedEvent(Node node) {
			this(node, null, -1);
		}

		public TreeChangedEvent(Node node, Parameter parameter, int index) {
			this.node = node;
			this.parameter = parameter;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		public Node getNode() {
			return node;
		}

		public Parameter getParameter() {
			return parameter;
		}

		public boolean isTreeChanged() {
			return parameter == null;
		}

		public boolean isNodeChanged() {
			return node != null;
		}

		public boolean isNodeParameterChanged() {
			return parameter != null;
		}

		public boolean isHeightChanged() {
			return parameter == node.heightParameter;
		}

		public boolean isRateChanged() {
			return parameter == node.rateParameter;
		}

		public boolean isTraitChanged(String name) {
			return parameter == node.traitParameters.get(name);
		}
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

	/**
	 * Sets the units that this tree is expressed in.
	 */
	public final void setUnits(Type units) {
		this.units = units;
	}

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
		return ((Node) node).getRate();
	}

	public Object getNodeAttribute(NodeRef node, String name) {
		throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
	}

	public Iterator getNodeAttributeNames(NodeRef node) {
		throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
	}

	public boolean hasNodeTraits() {
		return hasTraits;
	}

	public Map<String, Parameter> getTraitMap(NodeRef node) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		return ((Node) node).getTraitMap();
	}

	public double getNodeTrait(NodeRef node, String name) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		return ((Node) node).getTrait(name);
	}

	public Parameter getNodeTraitParameter(NodeRef node, String name) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		return ((Node) node).getTraitParameter(name);
	}

	public double[] getMultivariateNodeTrait(NodeRef node, String name) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		return ((Node) node).getMultivariateTrait(name);
	}

	public final void swapAllTraits(NodeRef node1, NodeRef node2) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		swapAllTraits((Node) node1, (Node) node2);
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
		return ((Node) node).getChildCount();
	}

	public final NodeRef getChild(NodeRef node, int i) {
		return ((Node) node).getChild(i);
	}

	public final NodeRef getParent(NodeRef node) {
		return ((Node) node).parent;
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

	public NodeRef[] getNodes() {
		return nodes;
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

	// *****************************************************************
	// Interface MutableTree
	// *****************************************************************

	/**
	 * Set a new node as root node.
	 */
	public final void setRoot(NodeRef newRoot) {

		if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

		root = (Node) newRoot;

		// We shouldn't need this because the addChild will already have fired appropriate events.
		// pushTreeChangedEvent();
	}

	public void addChild(NodeRef p, NodeRef c) {

		if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

		Node parent = (Node) p;
		Node child = (Node) c;
		if (parent.hasChild(child)) throw new IllegalArgumentException("Child already exists in parent");

		parent.addChild(child);
		pushTreeChangedEvent(parent);
	}

	public void removeChild(NodeRef p, NodeRef c) {

		if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

		Node parent = (Node) p;
		Node child = (Node) c;

		parent.removeChild(child);
	}

	private Node oldRoot;

	public void beginTreeEdit() {
		if (inEdit) throw new RuntimeException("Alreading in edit transaction mode!");

		oldRoot = root;

		inEdit = true;
	}

	public void endTreeEdit() throws MutableTree.InvalidTreeException {
		if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");

		inEdit = false;

		if (root != oldRoot) {
			swapParameterObjects(oldRoot, root);
		}

		for (Node node : nodes) {
			if (!node.heightParameter.isWithinBounds()) {
				throw new InvalidTreeException("height parameter out of bounds");
			}
		}

		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
			listenerHelper.fireModelChanged(this, treeChangedEvent);
		}
		treeChangedEvents.clear();
	}

	public void setNodeHeight(NodeRef n, double height) {
		((Node) n).setHeight(height);
	}


	public void setNodeRate(NodeRef n, double rate) {
		if (!hasRates) throw new IllegalArgumentException("Rate parameters have not been created");
		((Node) n).setRate(rate);

	}

	public void setNodeTrait(NodeRef n, String name, double value) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		((Node) n).setTrait(name, value);
	}

	public void setMultivariateTrait(NodeRef n, String name, double[] value) {
		if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
		((Node) n).setMultivariateTrait(name, value);
	}

	public void setBranchLength(NodeRef node, double length) {
		throw new UnsupportedOperationException("TreeModel cannot have branch lengths set");
	}

	public void setNodeAttribute(NodeRef node, String name, Object value) {
		throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
	}

	// *****************************************************************
	// Interface ModelComponent
	// *****************************************************************

	/**
	 * Store current state
	 */
	protected void storeState() {

		copyNodeStructure(storedNodes);
		storedRootNumber = root.getNumber();

	}

	/**
	 * Restore the stored state
	 */
	protected void restoreState() {

		Node[] tmp = storedNodes;
		storedNodes = nodes;
		nodes = tmp;

		root = nodes[storedRootNumber];
	}

	/**
	 * accept the stored state
	 */
	protected void acceptState() {
	} // nothing to do

	/**
	 * Copies the node connections from this TreeModel's nodes array to the
	 * destination array. Basically it connects up the nodes in destination
	 * in the same way as this TreeModel is set up. This method is package
	 * private.
	 */
	void copyNodeStructure(Node[] destination) {

		if (nodes.length != destination.length) {
			throw new IllegalArgumentException("Node arrays are of different lengths");
		}

		for (int i = 0, n = nodes.length; i < n; i++) {
			Node node0 = nodes[i];
			Node node1 = destination[i];

			// the parameter values are automatically stored and restored
			// just need to keep the links
			node1.heightParameter = node0.heightParameter;
			node1.rateParameter = node0.rateParameter;
			node1.traitParameters = node0.traitParameters;

			if (node0.parent != null) {
				node1.parent = storedNodes[node0.parent.getNumber()];
			} else {
				node1.parent = null;
			}

			if (node0.leftChild != null) {
				node1.leftChild = storedNodes[node0.leftChild.getNumber()];
			} else {
				node1.leftChild = null;
			}

			if (node0.rightChild != null) {
				node1.rightChild = storedNodes[node0.rightChild.getNumber()];
			} else {
				node1.rightChild = null;
			}
		}
	}

	/**
	 * @return the number of statistics of this component.
	 */
	public int getStatisticCount() {
		return 1;
	}

	/**
	 * @return the ith statistic of the component
	 */
	public Statistic getStatistic(int i) {
		if (i == 0) return root.heightParameter;
		throw new IllegalArgumentException();
	}

	public String getModelComponentName() {
		return TREE_MODEL;
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
		return Tree.Utils.newick(this);
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

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

	public Element createElement(Document document) {
		throw new RuntimeException("Not implemented yet");
	}

	// ***********************************************************************
	// Private methods
	// ***********************************************************************

	/**
	 * @return the node that this parameter is a member of
	 */
	public Node getNodeOfParameter(Parameter parameter) {

		if (parameter == null) throw new IllegalArgumentException("Parameter is null!");

		for (Node node : nodes) {
			if (node.heightParameter == parameter) {
				return node;
			}
			if (hasRates && node.rateParameter == parameter) {
				return node;
			}
		}
		for (Node node : nodes) {
			if (hasTraits && node.traitParameters.containsValue(parameter)) {
				return node;
			}
		}
		throw new RuntimeException("Parameter not found in any nodes:" + parameter.getId() + " " + parameter.hashCode());
		// assume it is a trait parameter and return null
//		return null;
	}

	/**
	 * Get the root height parameter. Is private because it can only be called by the XMLParser
	 */
	Parameter getRootHeightParameter() {

		return root.heightParameter;
	}

	/**
	 * @return the relevant node height parameter. Is private because it can only be called by the XMLParser
	 */
	Parameter createNodeHeightsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter("nodeHeights");

		for (int i = externalNodeCount; i < nodeCount; i++) {
			if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].heightParameter);
			}
		}

		if (leafNodes) {
			for (int i = 0; i < externalNodeCount; i++) {
				parameter.addParameter(nodes[i].heightParameter);
			}
		}

		return parameter;
	}

	Parameter getLeafHeightParameter(NodeRef node) {

		if (!isExternal(node)) {
			throw new RuntimeException("only root and leaves can be used with setNodeHeightParameter");
		}

		return nodes[node.getNumber()].heightParameter;
	}

	/**
	 * @return the relevant node rate parameter. Is private because it can only be called by the XMLParser
	 */
	Parameter createNodeRatesParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter("nodeRates");

		hasRates = true;

		for (int i = externalNodeCount; i < nodeCount; i++) {
			nodes[i].createRateParameter();
			if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].rateParameter);
			}
		}

		for (int i = 0; i < externalNodeCount; i++) {
			nodes[i].createRateParameter();
			if (leafNodes) {
				parameter.addParameter(nodes[i].rateParameter);
			}
		}

		return parameter;
	}

	/**
	 * Create a node traits parameter. Is private because it can only be called by the XMLParser
	 */
	public Parameter createNodeTraitsParameter(String name, int dim, double[] initialValues, boolean rootNode, boolean internalNodes, boolean leafNodes, boolean firesTreeEvents) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter("nodeTraits");

		hasTraits = true;

		for (int i = externalNodeCount; i < nodeCount; i++) {
			nodes[i].createTraitParameter(name, dim, initialValues, firesTreeEvents);
			if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].getTraitParameter(name));
			}
		}

		for (int i = 0; i < externalNodeCount; i++) {
			nodes[i].createTraitParameter(name, dim, initialValues, firesTreeEvents);
			if (leafNodes) {
				parameter.addParameter(nodes[i].getTraitParameter(name));
			}
		}

		return parameter;
	}

	private void swapAllTraits(Node n1, Node n2) {

		for (Map.Entry<String, Parameter> entry : n1.traitParameters.entrySet()) {
			Parameter p1 = n1.traitParameters.get(entry.getKey());
			Parameter p2 = n2.traitParameters.get(entry.getKey());
			final int dim = p1.getDimension();
			for (int i = 0; i < dim; i++) {
				double transfer = p1.getParameterValue(i);
				p1.setParameterValue(i, p2.getParameterValue(i));
				p2.setParameterValue(i, transfer);
			}

		}

	}

	/**
	 * This method swaps the parameter objects of the two nodes
	 * but maintains the values in each node.
	 * This method is used to ensure that root node of the tree
	 * always has the same parameter object.
	 */
	private void swapParameterObjects(Node n1, Node n2) {

		double height1 = n1.getHeight();
		double height2 = n2.getHeight();

		double rate1 = 1.0, rate2 = 1.0;

		if (hasRates) {
			rate1 = n1.getRate();
			rate2 = n2.getRate();
		}

		// swap all trait parameters

		if (hasTraits) {
			Map<String, Parameter> traits1 = new HashMap<String, Parameter>();
			Map<String, Parameter> traits2 = new HashMap<String, Parameter>();

			traits1.putAll(n1.traitParameters);
			traits2.putAll(n2.traitParameters);

			Map<String, Parameter> temp = n1.traitParameters;
			n1.traitParameters = n2.traitParameters;
			n2.traitParameters = temp;

			for (Map.Entry<String, Parameter> entry : traits1.entrySet()) {
				n1.traitParameters.get(entry.getKey()).setParameterValueQuietly(0, entry.getValue().getParameterValue(0));
			}
			for (Map.Entry<String, Parameter> entry : traits2.entrySet()) {
				n2.traitParameters.get(entry.getKey()).setParameterValueQuietly(0, entry.getValue().getParameterValue(0));
			}
		}

		Parameter temp = n1.heightParameter;
		n1.heightParameter = n2.heightParameter;
		n2.heightParameter = temp;

		if (hasRates) {
			temp = n1.rateParameter;
			n1.rateParameter = n2.rateParameter;
			n2.rateParameter = temp;
		}

		n1.heightParameter.setParameterValueQuietly(0, height1);
		n2.heightParameter.setParameterValueQuietly(0, height2);

		if (hasRates) {
			n1.rateParameter.setParameterValueQuietly(0, rate1);
			n2.rateParameter.setParameterValueQuietly(0, rate2);
		}

	}

	// **************************************************************
	// Private inner classes
	// **************************************************************

	public class Node implements NodeRef {

		public Node parent;
		public Node leftChild, rightChild;
		public int number;
		public Parameter heightParameter;
		public Parameter rateParameter = null;
		//public Parameter traitParameter = null;
		public Taxon taxon = null;

		Map<String, Parameter> traitParameters = new HashMap<String, Parameter>();

		public Node() {
			parent = null;
			leftChild = rightChild = null;
			heightParameter = null;
			number = 0;
			taxon = null;
		}

		/**
		 * constructor used to clone a node and all children
		 */
		public Node(Tree tree, NodeRef node) {
			parent = null;
			leftChild = rightChild = null;

			heightParameter = new Parameter.Default(tree.getNodeHeight(node));
			addParameter(heightParameter);

			number = node.getNumber();
			taxon = tree.getNodeTaxon(node);
            heightParameter.setId(""+number);
			for (int i = 0; i < tree.getChildCount(node); i++) {
				addChild(new Node(tree, tree.getChild(node, i)));
			}
		}

		public final void setupHeightBounds() {
			heightParameter.addBounds(new NodeHeightBounds(heightParameter));
		}

		public final void createRateParameter() {
			if (rateParameter == null) {
				rateParameter = new Parameter.Default(1.0);
				if (isRoot()) {
					rateParameter.setId("root.rate");
				} else if (isExternal()) {
					rateParameter.setId(getTaxonId(getNumber()) + ".rate");
				} else {
					rateParameter.setId("node" + getNumber() + ".rate");
				}
				rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
				addParameter(rateParameter);
			}
		}


		public final void createTraitParameter(String name, int dim, double[] initialValues, boolean firesTreeEvents) {

			if (!traitParameters.containsKey(name)) {

				Parameter trait = new Parameter.Default(dim);
				if (isRoot()) {
					trait.setId("root." + name);
				} else if (isExternal()) {
					trait.setId(getTaxonId(getNumber()) + "." + name);
				} else {
					trait.setId("node" + getNumber() + "." + name);
				}
				trait.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));

				if (initialValues != null && initialValues.length > 0) {
					for (int i = 0; i < dim; i++) {
						if (initialValues.length == dim) {
							trait.setParameterValue(i, initialValues[i]);
						} else {
							trait.setParameterValue(i, initialValues[0]);
						}
					}
				}

				traitParameters.put(name, trait);

				if (firesTreeEvents) {
					addParameter(trait);
				}
			}
		}

		public final double getHeight() {
			return heightParameter.getParameterValue(0);
		}

		public final double getRate() {
			return rateParameter.getParameterValue(0);
		}

		public final double getTrait(String name) {
			return traitParameters.get(name).getParameterValue(0);
		}

		public final double[] getMultivariateTrait(String name) {
			return traitParameters.get(name).getParameterValues();
		}

		public final Map<String, Parameter> getTraitMap() {
			return traitParameters;
		}

		public final void setHeight(double height) {
			heightParameter.setParameterValue(0, height);
		}

		public final void setRate(double rate) {
			//System.out.println("Rate set for parameter " + rateParameter.getParameterName());
			rateParameter.setParameterValue(0, rate);
		}

		public final void setTrait(String name, double trait) {
			//System.out.println("Trait set for parameter " + traitParameter.getParameterName());
			traitParameters.get(name).setParameterValue(0, trait);
		}

		public final void setMultivariateTrait(String name, double[] trait) {
			int dim = trait.length;
			for (int i = 0; i < dim; i++)
				traitParameters.get(name).setParameterValue(i, trait[i]);
		}

		public int getNumber() {
			return number;
		}

		/**
		 * Returns the number of children this node has.
		 */
		public final int getChildCount() {
			int n = 0;
			if (leftChild != null) n++;
			if (rightChild != null) n++;
			return n;
		}

		public Node getChild(int n) {
			if (n == 0) return leftChild;
			if (n == 1) return rightChild;
			throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
		}

		public boolean hasChild(Node node) {
			return (leftChild == node || rightChild == node);
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
				throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
			}
			node.parent = this;
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
			node.parent = null;
			return node;
		}

		/**
		 * remove child
		 *
		 * @param n number of child to be removed
		 */
		public Node removeChild(int n) {
			Node node;
			if (n == 0) {
				node = leftChild;
				leftChild = null;
			} else if (n == 1) {
				node = rightChild;
				rightChild = null;
			} else {
				throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
			}
			node.parent = null;
			return node;
		}

		public boolean hasChildren() {
			return (leftChild != null || rightChild != null);
		}

		public boolean isExternal() {
			return !hasChildren();
		}

		public boolean isRoot() {
			return (parent == null);
		}

		public String toString() {
			return "node " + number + ", height=" + getHeight() + (taxon != null ? ": " + taxon.getId() : "");
		}

		public Parameter getTraitParameter(String name) {
			return traitParameters.get(name);
		}
	}

	/**
	 * This class provides bounds for parameters that represent a node height
	 * in this tree model.
	 */
	private class NodeHeightBounds implements Bounds {

		public NodeHeightBounds(Parameter parameter) {
			nodeHeightParameter = parameter;
		}

		public double getUpperLimit(int i) {

			Node node = getNodeOfParameter(nodeHeightParameter);
			if (node.isRoot()) {
				return Double.POSITIVE_INFINITY;
			} else {
				return node.parent.getHeight();
			}
		}

		public double getLowerLimit(int i) {

			Node node = getNodeOfParameter(nodeHeightParameter);
			if (node.isExternal()) {
				return 0.0;
			} else {
				return Math.max(node.leftChild.getHeight(), node.rightChild.getHeight());
			}
		}

		public int getBoundsDimension() {
			return 1;
		}


		private Parameter nodeHeightParameter = null;
	}

	// ***********************************************************************
	// Private members
	// ***********************************************************************


	/**
	 * root node
	 */
	private Node root = null;
	private int storedRootNumber;

	/**
	 * list of internal nodes (including root)
	 */
	private Node[] nodes = null;
	private Node[] storedNodes = null;

	/**
	 * number of nodes (including root and tips)
	 */
	private int nodeCount;

	/**
	 * number of external nodes
	 */
	private int externalNodeCount;

	/**
	 * number of internal nodes (including root)
	 */
	private int internalNodeCount;

	/**
	 * holds the units of the trees branches.
	 */
	private Type units = Type.SUBSTITUTIONS;

	private boolean inEdit = false;

	private boolean hasRates = false;
	private boolean hasTraits = false;

}
