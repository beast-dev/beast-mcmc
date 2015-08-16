/*
 * UnrootedTree.java
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

package dr.evolution.tree;

public interface UnrootedTree {

	public int getNodeCount();

	public NodeRef getNode(int i);

	public int getNeighbourCount(NodeRef focal);

	public boolean isExternal(NodeRef focal);

	/**
	 * Get the branch length with the ith neighbour
	 */
	public double getBranchLength(NodeRef focal, int i);

	public void setBranchLength(NodeRef focal, int i, double edgeLength);

	/**
	 * Get the ith neighbouring node from the focal node
	 */
	public NodeRef getNeighbour(NodeRef focal, int i);

	/**
	 * Swap node A with node B across focal edge between node1 and node2
	 */
	public void NNI(NodeRef node1, int neighbour, int A, int B);


	public class Default implements UnrootedTree {


		public Default(Tree tree) {

			if (tree.getChildCount(tree.getRoot()) == 3) {

				nodeCount = tree.getNodeCount();

				nodeRef = new NodeRef[nodeCount];
				for (int i = 0; i < nodeCount; i++) {
					nodeRef[i] = new UnrootedNode(i);

					NodeRef treeNode = tree.getNode(i);
					if (!tree.isRoot(treeNode)) {
						NodeRef treeParent = tree.getParent(treeNode);
						addEdge(nodeRef[i], nodeRef[treeParent.getNumber()], tree.getNodeHeight(treeParent) - tree.getNodeHeight(treeNode));
					}
				}
			}

			if (tree.getChildCount(tree.getRoot()) == 2) {

				nodeCount = tree.getNodeCount() - 1;

				int count = 0;
				int[] mapping = new int[tree.getNodeCount()];
				for (int i = 0; i < tree.getNodeCount(); i++) {
					mapping[i] = count;
					if (!tree.isRoot(tree.getNode(i))) {
						count += 1;
					}
				}

				nodeRef = new NodeRef[nodeCount];
				for (int i = 0; i < tree.getNodeCount(); i++) {

					NodeRef treeNode = tree.getNode(i);
					if (!tree.isRoot(treeNode)) {

						nodeRef[mapping[i]] = new UnrootedNode(mapping[i]);

						NodeRef treeParent = tree.getParent(treeNode);

						if (!tree.isRoot(treeParent)) {
							addEdge(nodeRef[mapping[i]], nodeRef[mapping[treeParent.getNumber()]], tree.getNodeHeight(treeParent) - tree.getNodeHeight(treeNode));
						} else {

							NodeRef child1 = tree.getChild(tree.getRoot(), 0);
							NodeRef child2 = tree.getChild(tree.getRoot(), 1);

							double edgeLength =
								(tree.getNodeHeight(tree.getRoot()) - tree.getNodeHeight(child1)) +
								(tree.getNodeHeight(tree.getRoot()) - tree.getNodeHeight(child2));

							addEdge(nodeRef[mapping[child1.getNumber()]], nodeRef[mapping[child2.getNumber()]], edgeLength);
						}
					}
				}
			}
		}

		public Default(int tipCount) {

			nodeCount = tipCount * 2 - 1;

			edges = new double[nodeCount][nodeCount];

			for (int i = 0; i < nodeCount; i++) {
				for (int j = 0; j < nodeCount; j++) {
					edges[i][j] = -1.0;
				}
			}

			nodeRef = new NodeRef[nodeCount];
			for (int i = 0; i < nodeCount; i++) {
				nodeRef[i] = new UnrootedNode(i);
			}
		}

		public int getNodeCount() { return nodeCount; }
		public NodeRef getNode(int i) { return nodeRef[i];}

		/**
		 * @return the number of edges attached to the given node
		 */
		public int getNeighbourCount(NodeRef node) {
			int index = node.getNumber();
			int count = 0;
			for (int i = 0; i < nodeCount; i++) {
				if (edges[i][index] >= 0.0) {
					count += 1;
				}
			}
			return count;
		}

		/**
		 * @return true if this node has only one edge attached to it.
		 */
		public boolean isExternal(NodeRef node) {
			return getNeighbourCount(node) == 1;
		}

		public double getBranchLength(NodeRef node, int j) {
			int index = node.getNumber();
			int count = 0;
			for (int i = 0; i < nodeCount; i++) {
				if (edges[i][index] >= 0.0) {
					if (count == j) return edges[i][index];
					count += 1;
				}
			}
			throw new IllegalArgumentException();
		}

		public void setBranchLength(NodeRef node, int j, double edgeLength) {
			int index = node.getNumber();
			int count = 0;
			for (int i = 0; i < nodeCount; i++) {
				if (edges[i][index] >= 0.0) {
					if (count == j) {
						edges[i][index] = edgeLength;
						edges[index][i] = edgeLength;
					}
					count += 1;
				}
			}
		}

		/**
		 * @return the node at the end of the ith edge from the given node.
		 */
		public NodeRef getNeighbour(NodeRef node, int j) {
			int index = node.getNumber();
			int count = 0;
			for (int i = 0; i < nodeCount; i++) {
				if (edges[i][index] >= 0.0) {
					if (count == j) return nodeRef[i];
					count += 1;
				}
			}
			throw new IllegalArgumentException();
		}



		/**
		 * Swap node A with node B across focal edge between node1 and node2
		 */
		public void NNI(NodeRef node1, int neighbour, int A, int B) {

			if (neighbour == A) throw new IllegalArgumentException();

			NodeRef node2 = getNeighbour(node1, neighbour);

			NodeRef nA = getNeighbour(node1, A);
			NodeRef nB = getNeighbour(node2, B);

			addEdge(node1, nB, getBranchLength(node1, A));
			addEdge(node2, nA, getBranchLength(node2, B));
			removeEdge(node1, nA);
			removeEdge(node2, nB);
		}

		/*
		 * Should be called from root node.
		 * /
		private void traverse(NodeRef node, NodeRef from) {

			int edgeCount = getNeighbourCount(node);

			for (int i = 0; i < edgeCount; i++) {

				NodeRef neighbour = getNeighbour(node, i);

				if (from != neighbour) {
					traverse(neighbour, node);
				}
			}
		}*/

		/**
		 * Add an edge between the two nodes of the given length
		 */
		private void addEdge(NodeRef node, NodeRef node2, double edgeLength) {
			edges[node.getNumber()][node2.getNumber()] = edgeLength;
			edges[node2.getNumber()][node.getNumber()] = edgeLength;
		}

		/**
		 * Removes any edge between the two given nodes.
		 */
		private void removeEdge(NodeRef node, NodeRef node2) {
			edges[node.getNumber()][node2.getNumber()] = -1.0;
			edges[node2.getNumber()][node.getNumber()] = -1.0;
		}


		class UnrootedNode implements NodeRef {

			public UnrootedNode(int i) { number = i; }

			public int getNumber() { return number;}

                              public void setNumber(int n) { number = n; }

			int number;
		}

		private double[][] edges;
		private NodeRef[] nodeRef;
		private int nodeCount;
	}
}
