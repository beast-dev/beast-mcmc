/*
 * ARGModel.java
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

/*
 * ARGModel.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.arg;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.arg.likelihood.ARGLikelihood;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.inference.parallel.MPIServices;
import dr.math.MathUtils;
import dr.util.Attributable;
import dr.util.NumberFormatter;
import dr.xml.*;
import org.jdom.Document;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Logger;

/**
 * A model component for trees.
 *
 * @author Marc A. Suchard
 * @author Eric Bloomquist
 */
public class ARGModel extends AbstractModel implements MutableTree, Loggable {


    public static final String TREE_MODEL = "argTreeModel";
    public static final String ROOT_HEIGHT = TreeModelParser.ROOT_HEIGHT;
    public static final String LEAF_HEIGHT = "leafHeight";
    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_RATES = TreeModelParser.NODE_RATES;
    public static final String NODE_TRAITS = TreeModelParser.NODE_TRAITS;
    public static final String ROOT_NODE = "rootNode";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String LEAF_NODES = "leafNodes";
    public static final String TAXON = "taxon";
    public static final String GRAPH_ELEMENT = "graph";
    public static final String NODE_ELEMENT = "node";
    public static final String EDGE_ELEMENT = "edge";
    public static final String ID_ATTRIBUTE = XMLParser.ID;
    public static final String EDGE_FROM = "source";
    public static final String EDGE_TO = "target";
    public static final String TAXON_NAME = "taxonName";
    public static final String EDGE_LENGTH = "len";
    public static final String EDGE_PARTITIONS = "edgePartitions";
    public static final String IS_TIP = "isTip";
    public static final String IS_ROOT = "isRoot";
    public static final String LEFT_PARENT = "leftParent";
    public static final String RIGHT_PARENT = "rightParent";
    public static final String LEFT_CHILD = "leftChild";
    public static final String RIGHT_CHILD = "rightChild";
    public static final String NODE_HEIGHT = "nodeHeight";
    public static final String IS_REASSORTMENT = "true";
    public static final String NUM_PARTITIONS = "numberOfPartitions";
    public static final String GRAPH_SIZE = "size=\"6,6\"";
    public static final String DOT_EDGE_DEF = "edge[style=\"setlinewidth(2)\",arrowhead=none]";
    public static final String DOT_NODE_DEF = "node[shape=plaintext,width=auto,fontname=Helvitica,fontsize=10]";

    public static final int LEFT = 0;
    public static final int RIGHT = 1;

    public static final String PARTITION_TYPE = "partitionType";
    public static final String REASSORTMENT_PARTITION = "reassortment";
    public static final String RECOMBINATION_PARTITION = "recombination";
    public static final String PARTITION_DEFAULT_TYPE = REASSORTMENT_PARTITION;

    // ***********************************************************************
    // Private members
    // ***********************************************************************

    protected int storedRootNumber;
    protected int nodeCount;
    protected int storedNodeCount;
    protected int externalNodeCount;
    protected int internalNodeCount;
    protected int storedInternalNodeCount;
    protected boolean inEdit = false;

    protected Node root = null;
    public ArrayList<Node> nodes = null;
    public ArrayList<Node> storedNodes = null;
    protected Parameter[] addedParameters = null;
    protected Parameter[] removedParameters = null;
    protected Parameter addedPartitioningParameter = null;
    protected Parameter removedPartitioningParameter = null;
    protected CompoundParameter partitioningParameters;
    protected CompoundParameter storedInternalNodeHeights;
    protected CompoundParameter storedInternalAndRootNodeHeights;
    protected CompoundParameter storedNodeRates;
    protected Node[] addedNodes = null;
    protected Node[] removedNodes = null;

    //	private int units = SUBSTITUTIONS;
    private Type units;
    private boolean hasRates = false;
    private boolean hasTraits = false;
    private int nullCounter = 0;
    private int storedNullCounter;

    protected String partitionType = PARTITION_DEFAULT_TYPE;


    public ARGModel(ArrayList<Node> nodes, Node root, int numberPartitions,
                    int externalNodeCount) {
        super(TREE_MODEL);
        this.nodes = nodes;
        this.root = root;
        this.maxNumberOfPartitions = numberPartitions;
        this.externalNodeCount = externalNodeCount;
        if (nodes != null)
            this.nodeCount = nodes.size();
        this.internalNodeCount = this.nodeCount - externalNodeCount;
    }

    public ARGModel(Tree tree) {

        super(TREE_MODEL);
        // System.err.println("constructor for TreeModel");
        partitioningParameters = new CompoundParameter("partitioning");
        // initialize(tree);
        // }

        // protected void initialize(Tree tree) {
        // System.err.println("init for TreeModel");
        // System.exit(-1);
        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree);
        binaryTree.resolveTree();

        // clone the node structure (this will create the individual parameters
        Node node = new Node(binaryTree, binaryTree.getRoot());

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        // nodes = new Node[nodeCount];
        // storedNodes = new Node[nodeCount];
        nodes = new ArrayList<Node>(nodeCount);
        storedNodes = new ArrayList<Node>(nodeCount);
        for (int i = 0; i < nodeCount; i++) {
            nodes.add(null);
            storedNodes.add(null);
        }
        int i = 0;
        int j = externalNodeCount;

        root = node;

        // System.err.println("Going to do postOrder");

        do {
            node = (Node) TreeUtils.postorderSuccessor(this, node);

            if (node.isExternal()) {
                node.number = i;

                // nodes[i] = node;
                // storedNodes[i] = new Node();
                // storedNodes[i].taxon = node.taxon;
                // storedNodes[i].number = i;
                nodes.set(i, node);
                Node copy = new Node();
                copy.taxon = node.taxon;
                copy.number = i;
                storedNodes.set(i, copy);

                i++;
            } else {
                node.number = j;

                // nodes[j] = node;
                // storedNodes[j] = new Node();
                // storedNodes[j].number = j;
                nodes.set(j, node);
                Node copy = new Node();
                copy.number = j;
                storedNodes.set(j, copy);

                j++;
            }
        } while (node != root);

        // System.err.println("Succeed in post-order");

        // ARGTree t = new ARGTree(this,0);
        // System.err.println(Tree.Utils.uniqueNewick(t, t.getRoot()));
        // System.err.println(this.toGraphString());
        // System.exit(-1);

    }

    private double nextTime(int nTaxa, double pSize, double rRate) {
        double t = (double) nTaxa;

        double s = t * (t - 1 + rRate) / (2.0 * pSize);
        return MathUtils.nextExponential(s);
    }

    private boolean nextEventIsBifurcation(int nTaxa, double rRate) {
        double a = (double) (nTaxa - 1) / (nTaxa - 1 + rRate);
        if (MathUtils.nextDouble() < a) {
            return true;
        }
        return false;
    }

    private class SimulateSticks {
        public final Node mySon;
        public final boolean leftStick;

        public SimulateSticks(Node son, boolean left) {
            mySon = son;
            leftStick = left;
        }
    }

    public ARGModel(int ntaxa, double popSize, double rRate) {
        super("Simulator");

        ArrayList<SimulateSticks> currentStickList = new ArrayList<SimulateSticks>(50);
        ArrayList<Node> currentNodeList = new ArrayList<Node>(50);

        nodes = new ArrayList<Node>();

        int nodeNumber = 0;

        for (int i = 0; i < ntaxa; i++) {
            Node node = new Node();
            node.heightParameter = new Parameter.Default(0.0);
            nodes.add(node);
            currentNodeList.add(node);
            node.bifurcation = true;
            node.number = nodeNumber;
            nodeNumber++;

            node.taxon = new Taxon("" + nodeNumber);
            SimulateSticks stickGuy = new SimulateSticks(node, true);
            currentStickList.add(stickGuy);
        }

        double currentHeight = 0;

        while (currentStickList.size() > 1) {
            currentHeight = currentHeight + nextTime(currentStickList.size(), popSize, rRate);

            if (nextEventIsBifurcation(currentStickList.size(), rRate)) {
                SimulateSticks[] sticks = new SimulateSticks[2];

                for (int i = 0; i < 2; i++) {
                    int randomDraw = MathUtils.nextInt(currentStickList.size());
                    sticks[i] = currentStickList.get(randomDraw);
                    currentStickList.remove(randomDraw);
                }

                Node node = new Node();
                nodes.add(node);
                node.heightParameter = new Parameter.Default(currentHeight);
                node.number = nodeNumber;
                nodeNumber++;
                node.bifurcation = true;
                SimulateSticks stickGuy = new SimulateSticks(node, true);

                if (sticks[0].mySon == sticks[1].mySon) {
                    Node child = sticks[0].mySon;
                    node.leftChild = child;
                    node.rightChild = child;
                    child.rightParent = node;
                    child.leftParent = node;

                    currentNodeList.remove(currentNodeList.indexOf(child));
                    currentNodeList.add(node);
                    currentStickList.add(stickGuy);
                } else {
                    Node child1 = sticks[0].mySon;
                    Node child2 = sticks[1].mySon;

                    node.leftChild = child1;
                    node.rightChild = child2;

                    if (child1.bifurcation) {
                        child1.leftParent = node;
                        child1.rightParent = node;
                        currentNodeList.remove(currentNodeList.indexOf(child1));
                    } else {
                        if (sticks[0].leftStick) {
                            child1.leftParent = node;
                        } else {
                            child1.rightParent = node;
                        }
                        if (child1.leftParent != null && child1.rightParent != null) {
                            currentNodeList.remove(currentNodeList.indexOf(child1));
                        }
                    }

                    if (child2.bifurcation) {
                        child2.leftParent = node;
                        child2.rightParent = node;
                        currentNodeList.remove(currentNodeList.indexOf(child2));
                    } else {
                        if (sticks[1].leftStick) {
                            child2.leftParent = node;
                        } else {
                            child2.rightParent = node;
                        }
                        if (child2.leftParent != null && child2.rightParent != null) {
                            currentNodeList.remove(currentNodeList.indexOf(child2));
                        }
                    }
                    currentNodeList.add(node);
                    currentStickList.add(stickGuy);
                }
            } else {
                int randomDraw = MathUtils.nextInt(currentStickList.size());
                SimulateSticks stick = currentStickList.get(randomDraw);
                currentStickList.remove(randomDraw);

                Node node = new Node();
                nodes.add(node);
                node.heightParameter = new Parameter.Default(currentHeight);
                node.number = nodeNumber;
                nodeNumber++;
                node.bifurcation = false;
                SimulateSticks leftStickGuy = new SimulateSticks(node, true);
                SimulateSticks rightStickGuy = new SimulateSticks(node, false);

                Node child = stick.mySon;
                node.leftChild = child;
                node.rightChild = child;

                if (child.bifurcation) {
                    child.leftParent = node;
                    child.rightParent = node;
                    currentNodeList.remove(currentNodeList.indexOf(child));
                } else {
                    if (stick.leftStick) {
                        child.leftParent = node;
                    } else {
                        child.rightParent = node;
                    }
                    if (child.leftParent != null & child.rightParent != null) {
                        currentNodeList.remove(currentNodeList.indexOf(child));
                    }
                }

                currentNodeList.add(node);
                currentStickList.add(leftStickGuy);
                currentStickList.add(rightStickGuy);
            }

        }

        root = currentNodeList.get(0);

        //nodeNumber--;

        //	nodes = new Node[nodeNumber];


    }


    public int possibleInternalNodePermuations() {

        int max = getInternalNodeCount();
        ArrayList<Double> remainingHeights = new ArrayList<Double>(max - 1);

        for (Node node : nodes) {
            if (!node.isExternal() && !node.isRoot())
                remainingHeights.add(node.getHeight());
        }

        int firstNode = 0;
        while (nodes.get(firstNode).isExternal() || nodes.get(firstNode).isRoot())
            firstNode++;

        int result = possibleInternalNodePermutations(firstNode, remainingHeights);

        // Restore heights
        int i = 0;
        for (Node node : nodes) {
            if (!node.isExternal() && !node.isRoot()) {
                node.setHeight(remainingHeights.get(i++));
            }
        }

        return factorial(max - 1) - result;
    }

    private int factorial(int x) {
        int result = 1;
        for (int i = 2; i <= x; i++)
            result *= i;
        return result;
    }


    private int possibleInternalNodePermutations(int nodeNumber, ArrayList<Double> remainingHeights) {

        int total = 0;

        if (remainingHeights.size() == 0) {
            return 0;
        }

        int newNodeNumber = nodeNumber + 1;
        if (remainingHeights.size() > 1) {
            while (nodes.get(newNodeNumber).isExternal() || nodes.get(newNodeNumber).isRoot())
                newNodeNumber++;
        }

        Node nr = nodes.get(nodeNumber);

        for (double height : remainingHeights) {

            if (height < getNodeHeight(getParent(nr, 0)) &&
                    height < getNodeHeight(getParent(nr, 1))) {
                setNodeHeight(nr, height);

                ArrayList<Double> copy0 = deepCopy(remainingHeights);
                if (!copy0.contains(height))
                    System.err.println("where did i go?");
                copy0.remove(height);

                total += possibleInternalNodePermutations(newNodeNumber, copy0);

            } else {
                // The remaining permutations will not work
                total += factorial(remainingHeights.size() - 1);
            }

        }

        return total;
    }


    private ArrayList<Double> deepCopy(ArrayList<Double> in) {
        ArrayList<Double> out = new ArrayList<Double>();
        for (double d : in) {
            out.add(d);
        }
        return out;
    }

    private static boolean containsLessThan(int[] a, int b) {
        for (int c : a) {
            if (c < b)
                return true;
        }
        return false;
    }

    public static void main(String[] args) {
        ARGModel arg = new ARGModel(8, 20.0, 0.5);
        System.out.println(arg.toARGSummary());
        System.out.println(arg.getReassortmentNodeCount());
        System.out.println(arg.toExtendedNewick());
//		BufferedWriter out = null;
//
//		try {
//			out = new BufferedWriter(new FileWriter("rejections.txt"));
//		} catch (Exception e) {
//			System.exit(-1);
//		}
//
//		int rejections = 0;
//
//		MathUtils.setSeed(97695);
//
//		ArrayList<String> trees = new ArrayList<String>(10000);
//
//		while (rejections < 100) {
//			ARGModel arg = new ARGModel(3, 6.0, 1.0);
//			String s = "";
//			if (arg.getReassortmentNodeCount() < 2) {
//				s = arg.toStrippedNewick();
//			}
//
//			if (!s.equals("") && !trees.contains(s)) {
//				trees.add(s);
//				System.out.println("Total Current Size = " + trees.size() + " Rejections= " + rejections);
//				rejections = 0;
//			} else {
//				rejections++;
//			}
//		}
//
//		System.out.println("\n************************************");
//		System.out.println("Simulating trees");
//
//		Collections.sort(trees);
//		int[] freq = new int[trees.size()];
//		int[] mcmcFreq = new int[trees.size()];
//		rejections = 0;
//
//		while (containsLessThan(freq, 10000)) {
//			ARGModel arg = new ARGModel(3, 6.0, 1.0);
//			String s = null;
//			if (arg.getReassortmentNodeCount() < 2) {
//				s = arg.toStrippedNewick();
//			}
//
//			if (s != null) {
//				freq[Collections.binarySearch(trees, s)]++;
//			}
//			rejections++;
//			if (rejections % 1000000 == 0) {
//				System.out.println(rejections);
//			}
//		}
//
////		System.out.println("\n************************************");
////		System.out.println("Analyzing MCMC results");
////
////		BufferedReader read = null;
////		try{
////			read = new BufferedReader( new FileReader("prior.args"));
////			String s = read.readLine();
////			s = read.readLine();
////
////			while(s != null){
////				mcmcFreq[Collections.binarySearch(trees, s)]++;
////				s = read.readLine();
////			}
////
////		}catch(Exception e){
////			System.exit(-1);
////		}
//
//		System.out.println("\n************************************");
//		System.out.println("Printing Results");
//
//
//		try {
//			out = new BufferedWriter(new FileWriter("coalescent.sim"));
//
//			rejections = 0;
//			for (String s : trees) {
//				if (mcmcFreq[rejections] > -1) {
//					out.write(s + " " + freq[rejections] + " " + mcmcFreq[rejections] + " \n");
//				}
//				rejections++;
//			}
//			out.flush();
//		} catch (Exception IOException) {
//			System.exit(-1);
//		}

    }

    /**
     * Packs and sends ARG state, including connectedness, heightparameters and
     * partitioning parameters.
     *
     * @param toRank
     */

    @Override
    public void sendState(int toRank) {
        sendStateNoParameters(toRank);
        int cnt = 0;
        for (Node node : nodes) {
            node.number = cnt;
            cnt++;
        }
        final int size = nodes.size();
        MPIServices.sendInt(size, toRank);
        int[] intMsg = new int[size * 7];
        double[] doubleMsg = new double[size];
        int indexNode = 0;
        int indexHeight = 0;
        int indexPartition = 0;
        ArrayList<Parameter> partList = new ArrayList<Parameter>();
        for (Node node : nodes) {
            intMsg[indexNode++] = node.number;

            if (node.leftParent != null)
                intMsg[indexNode++] = node.leftParent.number;
            else
                intMsg[indexNode++] = -1;
            if (node.rightParent != null)
                intMsg[indexNode++] = node.rightParent.number;
            else
                intMsg[indexNode++] = -1;
            if (node.leftChild != null)
                intMsg[indexNode++] = node.leftChild.number;
            else
                intMsg[indexNode++] = -1;
            if (node.rightChild != null)
                intMsg[indexNode++] = node.rightChild.number;
            else
                intMsg[indexNode++] = -1;

            if (node.partitioning != null) {
                intMsg[indexNode++] = indexPartition++;
                partList.add(node.partitioning);
            } else {
                intMsg[indexNode++] = -1;
            }

            if (node.bifurcation)
                intMsg[indexNode++] = 1;
            else
                intMsg[indexNode++] = 0;

            doubleMsg[indexHeight++] = node.heightParameter
                    .getParameterValue(0);
        }
        MPIServices.sendIntArray(intMsg, toRank);
        MPIServices.sendDoubleArray(doubleMsg, toRank);
        MPIServices.sendInt(partList.size(), toRank);
        for (Parameter partition : partList) {
            // System.err.println("Sending a partition.");
            double[] values = partition.getParameterValues();
            // System.err.println("length = "+values.length);
            MPIServices.sendDoubleArray(partition.getParameterValues(), toRank);
        }
        // partitioningParameters.sendState(toRank);
        MPIServices.sendInt(((Node) getRoot()).number, toRank);
    }

    @Override
    public void receiveState(int fromRank) {
        receiveStateNoParameters(fromRank);
        final int newNodeCount = MPIServices.receiveInt(fromRank);
        // while (newNodeCount < nodes.size())
        // nodes.remove(0);
        int[] intMsg = MPIServices.receiveIntArray(fromRank, newNodeCount * 7);
        double[] doubleMsg = MPIServices.receiveDoubleArray(fromRank,
                newNodeCount);
        int partitionLength = MPIServices.receiveInt(fromRank);
        // System.err.println("Attemping to receive "+partitionLength+"
        // partitions.");
        final int length = getNumberOfPartitions();
        // System.err.println("Expected length = "+length);
        while (partitionLength > partitioningParameters.getParameterCount()) {
            Parameter newPartition = new Parameter.Default(length);
            partitioningParameters.addParameter(newPartition);
        }

        for (int i = 0; i < partitionLength; i++) {
            double[] values = MPIServices.receiveDoubleArray(fromRank, length);
            // System.err.println("Received.");
            Parameter param = partitioningParameters.getParameter(i);
            // System.err.println("null? "+ (param == null ? "Yes" : "No"));
            for (int j = 0; j < length; j++) {
                // System.err.println("setting value #"+j+ " = "+values[j]);

                param.setParameterValueQuietly(j, values[j]);
            }
            // System.err.println("Values set");
        }
        // System.err.println("Done with partition receive.");

        int root = MPIServices.receiveInt(fromRank);
        // System.err.println("Start reconstructing ARG");

        beginTreeEdit();
        while (newNodeCount > nodes.size()) {
            Node newNode = new Node();
            newNode.heightParameter = new Parameter.Default(0.0);
            nodes.add(newNode);
        }
        // System.err.println("extra height added");
        int nodeInt;
        int indexNode = 0;
        int indexHeight = 0;
        int indexPartition = 0;
        for (int i = 0; i < newNodeCount; i++) {
            Node node = nodes.get(i);
            node.number = intMsg[indexNode++];

            nodeInt = intMsg[indexNode++];
            if (nodeInt != -1)
                node.leftParent = nodes.get(nodeInt);
            else
                node.leftParent = null;
            nodeInt = intMsg[indexNode++];
            if (nodeInt != -1)
                node.rightParent = nodes.get(nodeInt);
            else
                node.rightParent = null;
            nodeInt = intMsg[indexNode++];
            if (nodeInt != -1)
                node.leftChild = nodes.get(nodeInt);
            else
                node.leftChild = null;
            nodeInt = intMsg[indexNode++];
            if (nodeInt != -1)
                node.rightChild = nodes.get(nodeInt);
            else
                node.rightChild = null;
            // Parameter heightParam =
            node.heightParameter.setParameterValueQuietly(0,
                    doubleMsg[indexHeight++]);
            int whichPartitionParameter = intMsg[indexNode++];
            if (whichPartitionParameter != -1) {
                // System.err.println("Setting partition para");

                node.partitioning = partitioningParameters
                        .getParameter(whichPartitionParameter);
                // System.err.println("Done setting param");
            }

            if (intMsg[indexNode++] == 1)
                node.bifurcation = true;
            else
                node.bifurcation = false;

        }
        // System.err.println("Recovered all nodes");
        setRoot(nodes.get(root));
        // try {
        endTreeEditFast();
        // todo fire an ARG changed event???

        // } catch (InvalidTreeException e) {
        // throw new RuntimeException("Unable to unpack ARG correctly");
        // e.printStackTrace(); //To change body of catch statement use File |
        // Settings | File Templates.
        // }

    }

    public boolean isAncestral() {

        //zero everything first.
        for (int i = 0; i < getNodeCount(); i++) {
            Node x = (Node) getNode(i);

            x.fullAncestralMaterial = false;
            x.hasSomeAncestralMaterial = false;

            if (x.ancestralMaterial == null) {
                x.ancestralMaterial = new boolean[getNumberOfPartitions()];
            }

            for (int j = 0; j < x.ancestralMaterial.length; j++) {
                x.ancestralMaterial[j] = false;
            }

        }

        //post order up with the external nodes

        for (int i = 0; i < getExternalNodeCount(); i++) {
            Node currentNode = (Node) getExternalNode(i);

            currentNode.fullAncestralMaterial = true;
            currentNode.hasSomeAncestralMaterial = true;

            for (int j = 0; j < currentNode.ancestralMaterial.length; j++) {
                currentNode.ancestralMaterial[j] = true;
            }

            currentNode.leftParent.setAncestralMaterial(currentNode.ancestralMaterial);
        }


        //check that everything has some ancestral stuff.
        for (int i = 0, n = getNodeCount(); i < n; i++) {
            Node currentNode = (Node) getNode(i);

            if (!currentNode.hasSomeAncestralMaterial) {
                return false;
            }
        }


        return true;
    }


    public CompoundParameter getPartitioningParameters() {
        return partitioningParameters;
    }

    public void setupHeightBounds() {
        for (Node node : nodes) {
            node.setupHeightBounds();
        }
    }

    // public static final String

    private String getNameOfNode(Node node) {
        if (node.taxon == null)
            return "n" + Integer.toString(node.number);
        else
            return node.taxon.getId();
        // return Integer.toString(node.number);
    }

    public static final int MAX_LABEL_COUNT = 10;

    private Element makeEdge(Node from, Node to) {
        Element edgeElement = new Element(EDGE_ELEMENT);
        edgeElement.setAttribute(EDGE_FROM, getNameOfNode(from));
        edgeElement.setAttribute(EDGE_TO, getNameOfNode(to));
        edgeElement.setAttribute(EDGE_LENGTH, Double
                .toString(getNodeHeight(from) - getNodeHeight(to)));
        if (to.isReassortment()) {
            double[] bits = to.partitioning.getParameterValues();
            int length = bits.length;
            StringBuilder sb = new StringBuilder();
            boolean isLeft = (from == to.leftParent);
            int countLabels = 0;
            for (int i = 0; i < length; i++) {
                if ((isLeft && bits[i] == 0) || (!isLeft && bits[i] == 1)) {
                    sb.append(i);
                    sb.append(" ");
                    countLabels++;
                }

            }
            if (countLabels < MAX_LABEL_COUNT)
                edgeElement.setAttribute(EDGE_PARTITIONS, sb.toString().trim());
        }
        return edgeElement;
    }

    private Element makeNode(Node node) {
        Element nodeElement = new Element(NODE_ELEMENT);
        nodeElement.setAttribute(ID_ATTRIBUTE, getNameOfNode(node));
        if (node.taxon != null) {
            nodeElement.setAttribute(IS_TIP, "true");
            nodeElement.setAttribute(TAXON_NAME, node.taxon.getId());
            // nodeElement.setAttribute("style","filled");
            // nodeElement.setAttribute("fillcolor","blue");
        }
        if (node.isRoot()) {
            nodeElement.setAttribute(IS_ROOT, "true");
        }
        if (node.isReassortment()) {
            nodeElement.setAttribute(IS_REASSORTMENT, "true");
        }
        nodeElement
                .setAttribute(NODE_HEIGHT, Double.toString(node.getHeight()));
        return nodeElement;
    }

    private Element makeNodeFullInfo(Node node) {
        Element nodeElement = new Element(NODE_ELEMENT);
        nodeElement.setAttribute(ID_ATTRIBUTE, getNameOfNode(node));

        if (node.isRoot()) {
            nodeElement.setAttribute(IS_ROOT, "true");
        } else {
            nodeElement.setAttribute(LEFT_PARENT,
                    getNameOfNode(node.leftParent));
            nodeElement.setAttribute(RIGHT_PARENT,
                    getNameOfNode(node.rightParent));
        }

        if (node.taxon != null) {
            nodeElement.setAttribute(IS_TIP, "true");
            nodeElement.setAttribute(TAXON_NAME, node.taxon.getId());
        } else {
            nodeElement.setAttribute(LEFT_CHILD, getNameOfNode(node.leftChild));
            nodeElement.setAttribute(RIGHT_CHILD,
                    getNameOfNode(node.rightChild));
        }
        nodeElement
                .setAttribute(NODE_HEIGHT, Double.toString(node.getHeight()));

        return nodeElement;
    }

    public ARGModel fromXML(Element rootElement) {

        int numPartitions = Integer.parseInt(rootElement
                .getAttributeValue(NUM_PARTITIONS));
        int external = 0;

        // count # of nodeElements
        List<Element> nodeList = rootElement.getChildren(NODE_ELEMENT);
        ArrayList<Node> nodes = new ArrayList<Node>();
        Node rootNode = null;
        // System.err.println("node # = "+nodeList.size());
        // System.exit(-1);
        for (Element nodeElement : nodeList) {
            Node node = new Node();
            nodes.add(node);
            String isRoot = nodeElement.getAttributeValue(IS_ROOT);
            if (isRoot != null && isRoot.compareTo("true") == 0)
                rootNode = node;
            String isReassortment = nodeElement
                    .getAttributeValue(IS_REASSORTMENT);
            if (isReassortment != null && isReassortment.compareTo("true") == 0)
                node.bifurcation = false;
            else
                node.bifurcation = true;
            double height = Double.parseDouble(nodeElement
                    .getAttributeValue(NODE_HEIGHT));

            node.heightParameter = new Parameter.Default(height);
            node.setHeight(height);

            String isTip = nodeElement.getAttributeValue(IS_TIP);
            if (isTip != null && isTip.compareTo("true") == 0) {
                external++;
                String taxonName = nodeElement.getAttributeValue(TAXON_NAME);
                node.taxon = new Taxon(taxonName); // todo reuse taxonList
            }
        }

        List<Element> edgeList = rootElement.getChildren(EDGE_ELEMENT);
        for (Element edgeElement : edgeList) {
            int target = Integer.parseInt(edgeElement
                    .getAttributeValue(EDGE_TO));
            int source = Integer.parseInt(edgeElement
                    .getAttributeValue(EDGE_FROM));
            Node targetNode = nodes.get(target);
            Node sourceNode = nodes.get(source);

            if (targetNode.isBifurcation()) {
                targetNode.leftParent = targetNode.rightParent = sourceNode;
            } else {
                if (targetNode.leftParent == null)
                    targetNode.leftParent = sourceNode;
                else {
                    targetNode.rightParent = sourceNode;
                    String partitionInfo = edgeElement
                            .getAttributeValue(EDGE_PARTITIONS);
                    Parameter partitioning = null;
                    // if (partitionInfo != null && sourceNode.leftChild !=
                    // null) {
                    partitioning = new Parameter.Default(numPartitions, 0.0);
                    StringTokenizer st = new StringTokenizer(partitionInfo);
                    while (st.hasMoreTokens()) {
                        int which = Integer.parseInt(st.nextToken());
                        partitioning.setParameterValueQuietly(which, 1.0);
                    }
                    // }
                    targetNode.partitioning = partitioning;

                }
            }

            if (sourceNode.leftChild == null)
                sourceNode.leftChild = targetNode;
                // sourceNode.
                // }
            else {
                sourceNode.rightChild = targetNode;
                // sourceNode.partitioning = partitioning;
            }
            // todo parse partition info
        }

        return new ARGModel(nodes, rootNode, numPartitions, external);
    }

    public Element toXML() {

        // toGraphStringCompressed(true);

        int cnt = 0;
        for (Node node : nodes)
            node.number = cnt++;

        Element graphElement = new Element(GRAPH_ELEMENT);
        graphElement.setAttribute("edgedefault", "directed");
        graphElement.setAttribute(NUM_PARTITIONS, Integer
                .toString(getNumberOfPartitions()));

        for (Node node : nodes) {

            graphElement.addContent(makeNode(node));
            // Add edge to left parent if not root
            if (node.leftParent != null) {
                graphElement.addContent(makeEdge(node.leftParent, node));

            }
            // Add edge to right parent if reassortment
            if (node.rightParent != null && node.isReassortment()) {
                graphElement.addContent(makeEdge(node.rightParent, node));
            }

        }

        // System.err.println("start = "+nodes.size());
        // ARGModel test = fromXML(graphElement);
        // System.err.println(test.toGraphString());

        // System.err.println("old 0:"+getNewick(0));
        // System.err.println("old 1:"+getNewick(1));
        // System.err.println("old 2:"+getNewick(2));
        // System.err.println("new 0:"+test.getNewick(0));

        return graphElement;
    }

    public String toStrippedNewick() {
        String s = root.toExtendedNewick() + ";";
//		s = s.replaceAll("[0-9a-zA-Z]", "");
        s = s.replaceAll("[^(),<>;]", "");
        return s;

    }

    public String toExtendedNewick() {
        // StringBuffer sb = new StringBuffer();
        return root.toExtendedNewick() + ";";
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent() {
        pushTreeChangedEvent(new ARGTreeChangedEvent());
    }

    public void pushTreeSizeChangedEvent() {
        throw new RuntimeException("No longer supported; use updated operators");
    }

    public void pushTreeSizeIncreasedEvent() {
        pushTreeChangedEvent(new ARGTreeChangedEvent(+1));
    }

    public void pushTreeSizeDecreasedEvent() {
        pushTreeChangedEvent(new ARGTreeChangedEvent(-1));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(NodeRef nodeRef) {
        pushTreeChangedEvent(new ARGTreeChangedEvent((Node) nodeRef));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(Node node, Parameter parameter, int index) {
        pushTreeChangedEvent(new ARGTreeChangedEvent(node, parameter, index));
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
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        Node node = getNodeOfParameter((Parameter) variable);
        pushTreeChangedEvent(node, (Parameter) variable, index);
    }

    private ArrayList<ARGLikelihood> likelihoodCalculators;

    private int maxNumberOfPartitions;

    public int getNumberOfPartitions() {
        return maxNumberOfPartitions;
    }

    public int addLikelihoodCalculator(ARGLikelihood calc) {
        // int len = 0;
        if (likelihoodCalculators == null) {
            likelihoodCalculators = new ArrayList<ARGLikelihood>();
        }
        likelihoodCalculators.add(calc);
        int len = likelihoodCalculators.size() - 1;
        maxNumberOfPartitions = likelihoodCalculators.size();
        System.err.println("Add calculator for partition #" + len);
        setPartitionRecursively(getRoot(), len);
        return len;
    }

    public int getMaxPartitionNumber() {
        return maxNumberOfPartitions;
    }

    protected final List<TreeChangedEvent> treeChangedEvents = new ArrayList<TreeChangedEvent>();

    public class ARGTreeChangedEvent implements TreeChangedEvent {

        final Node node;

        final Parameter parameter;

        final int index;

        int size = 0;

        public ARGTreeChangedEvent() {
            this(null, null, -1);
        }

//        public TreeChangedEvent(ARGModel arg) {
//            this(null, null, -1);
//            size = true;
//        }

        public ARGTreeChangedEvent(Node node) {
            this(node, null, -1);
        }

        public ARGTreeChangedEvent(Node node, Parameter parameter, int index) {
            this.node = node;
            this.parameter = parameter;
            this.index = index;
        }

        public ARGTreeChangedEvent(int sizeChanged) {
            this(null, null, -1);
            size = sizeChanged;
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

        public int getSize() {
            return size;
        }

        public boolean isSizeChanged() {
            return !(size == 0);
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

        public boolean isTraitChanged() {
            return parameter == node.traitParameter;
        }
    }

    // *****
    // Interface Loggable
    // *****

    public LogColumn[] getColumns() {
        int numColumns = 3;
        // numColumns += this.getMaxPartitionNumber();
        // LogColumn[] logColumns = new LogColumn[numColumns +
        // getMaxPartitionNumber()];
        LogColumn[] logColumns = new LogColumn[3];
        //logColumns[0] = new IsReassortmentColumn("isReassortment");
        logColumns[0] = new CountReassortmentColumn("numberReassortments");
        logColumns[1] = new ExtremeNodeHeightColumn("maxNodeHeight") {
            double getStartValue() {
                return 0;
            }

            double compare(double currentValue, double newValue) {
                if (newValue > currentValue)
                    return newValue;
                return currentValue;
            }
        };
        logColumns[2] = new ExtremeNodeHeightColumn("minNodeHeight") {
            double getStartValue() {
                return 0;
            }

            double compare(double currentValue, double newValue) {
                if (newValue == 0) {
                    return currentValue;
                } else if (currentValue == 0) {
                    return newValue;
                }
                if (newValue < currentValue)
                    return newValue;
                return currentValue;
            }
        };
        /*logColumns[4] = new NumberColumn("Left Node"){
              public double getDoubleValue() {
                  Node a = (Node) getRoot();
                  return a.leftChild.getHeight();
              }
          };
          logColumns[5] = new NumberColumn("Right Node"){
              public double getDoubleValue() {
                  Node a = (Node) getRoot();
                  return a.rightChild.getHeight();
              }
          };
          logColumns[6] = new NumberColumn("Reassort Height"){
              public double getDoubleValue() {
                  double b = 0;
                  for(Node a : nodes){
                      if(a.isReassortment()){
                          b = a.getHeight();
                          break;
                      }
                  }
                  return b;
              }
          };*/

        // logColumns[2] = new IsRootTooHighColumn("isRootTooHigh");
        // for (int i = 0; i < getMaxPartitionNumber(); i++)
        // logColumns[4 + i] = new ArgTreeHeightColumn("argTreeHeight", this,
        // i);
        return logColumns;
    }

    private abstract class ExtremeNodeHeightColumn extends NumberColumn {

        public ExtremeNodeHeightColumn(String label) {
            super(label);
        }

        abstract double compare(double currentValue, double newValue);

        /*
           * { if (newValue > currentValue) return newValue; return currentValue; }
           */

        abstract double getStartValue();

        public double getDoubleValue() {
            double criticalValue = getStartValue();
            for (Node node : nodes) {
                double nodeHeight = node.heightParameter.getParameterValue(0);
                if (nodeHeight > 0 && !isRoot(node)) {
                    criticalValue = compare(criticalValue, nodeHeight);
                }
            }
            return criticalValue;
        }
    }

    private class ArgTreeHeightColumn extends NumberColumn {

        private final int partition;

        private final ARGModel argModel;

        public ArgTreeHeightColumn(String label, ARGModel argModel,
                                   int partition) {
            super(label + partition);
            this.argModel = argModel;
            this.partition = partition;
        }

        public double getDoubleValue() {
            ARGTree argTree = new ARGTree(argModel, partition);
            return argTree.getNodeHeight(argTree.getRoot());
            // return (new ARGTree(
        }
    }

    private class IsReassortmentColumn extends NumberColumn {

        public IsReassortmentColumn(String label) {
            super(label); // To change body of overridden methods use File |
            // Settings | File Templates.
        }

        public double getDoubleValue() {
            return getReassortmentNodeCount() == 0 ? 0 : 1;
        }
    }

    private class IsRootTooHighColumn extends NumberColumn {

        public IsRootTooHighColumn(String label) {
            super(label); // To change body of overridden methods use File |
            // Settings | File Templates.
        }

        public double getDoubleValue() {
            return isBifurcationDoublyLinked(getRoot()) ? 1 : 0;
        }
    }

    private class CountReassortmentColumn extends NumberColumn {

        public CountReassortmentColumn(String label) {
            super(label); // To change body of overridden methods use File |
            // Settings | File Templates.
        }

        public double getDoubleValue() {
            return getReassortmentNodeCount();
        }
    }

    public void setPartitionType(String partitionType) {
        this.partitionType = partitionType;
    }

    public String getPartitionType() {
        return partitionType;
    }

    public boolean isRecombinationPartitionType() {
        if (partitionType.equals(RECOMBINATION_PARTITION)) {
            return true;
        }
        return false;
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
        return nodes.size();
    }

    public final boolean hasNodeHeights() {
        return true;
    }

    public NodeRef getMirrorNode(NodeRef node) {

//        for(Node argNode: nodes) {
//            if (((Node)node).mirrorNode == argNode)
//                System.err.println("Found (in getMirrorNode)");
//        }

        return ((Node) node).mirrorNode;
    }

    public final double getNodeHeight(NodeRef node) {

        // System.err.println(Tree.Utils.uniqueNewick(this, node));
        // ((Node)node))

        return ((Node) node).getHeight();
    }

    public final double getMinParentNodeHeight(NodeRef nr) {
        Node node = (Node) nr;
        return Math.min(node.leftParent.getHeight(), node.rightParent
                .getHeight());
    }

    public final double getNodeHeightUpper(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getUpperLimit(0);
    }

    public final double getNodeHeightLower(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getLowerLimit(0);
    }

    /**
     * @param nodeRef
     * @return the rate parameter associated with this node.
     */
    public final double getNodeRate(NodeRef nodeRef, int partition) {
        if (!hasRates) {
            return 1.0;
        }


        return 0.0;
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException(
                "ARGModel does not use NodeAttributes");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException(
                "ARGModel does not use NodeAttributes");
    }

    public double getNodeTrait(NodeRef node) {
        if (!hasTraits)
            throw new IllegalArgumentException(
                    "Trait parameters have not been created");
        return ((Node) node).getTrait();
    }

    public final Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public final boolean isExternal(NodeRef node) {
        return ((Node) node).isExternal();
    }

    public final boolean isInternal(NodeRef node) {
        return !this.isExternal(node);
    }

    public final boolean isRoot(NodeRef node) {
        return (node == root);
    }

    public final boolean isBifurcation(NodeRef node) {
        return ((Node) node).isBifurcation();
    }

    public final boolean isBifurcationDoublyLinked(NodeRef node) {
        return ((Node) node).isBifurcationDoublyLinked();
    }

    public final boolean isReassortment(NodeRef node) {
        return ((Node) node).isReassortment();
    }

    public final int countReassortmentNodes(NodeRef nr) {
        Node node = (Node) nr;
        int count = node.countReassortmentChild(this);
        // int count = 0;
        return (count / 2);
    }

    public final int getChildCount(NodeRef node) {
        return ((Node) node).getChildCount();
    }


    public final NodeRef getOtherChild(NodeRef parent, NodeRef wrongChild) {
        Node p = (Node) parent;
        Node c = (Node) wrongChild;

        if (p.leftChild == c) {
            return p.rightChild;
        }
        return p.leftChild;
    }

    public final NodeRef getBrother(NodeRef node) {
        Node n = (Node) node;

        if (n.isReassortment()) {
            return node;
        }
        Node p = n.leftParent;

        if (p.leftChild == n) {
            return p.rightChild;
        }

        return p.leftChild;
    }

    /**
     * If i = 0, the left child is returned, else if i = 1, the right child is
     * returned.
     *
     * @return The child of the entered node.
     */
    public final NodeRef getChild(NodeRef node, int i) {
        return ((Node) node).getChild(i);
    }

    public final NodeRef getChild(NodeRef node, int i, int partition) {
        return ((Node) node).getChild(i, partition);
    }

    // public final NodeRef getParent(NodeRef node) { return
    // ((Node)node).parent; }
    public final NodeRef getParent(NodeRef node) {
        Node left = ((Node) node).leftParent;
        Node right = ((Node) node).rightParent;
        if (left == right)
            return left;
        else
            throw new IllegalArgumentException(
                    "No single parent for reassorted node");
    }

    /**
     * @param node The child noderef
     * @param i    i = 0 (left parent) i = 1 (right parent)
     * @return The corresponding parent noderef
     */
    public final NodeRef getParent(NodeRef node, int i) {
        if (i == 0)
            return ((Node) node).leftParent;
        if (i == 1)
            return ((Node) node).rightParent;
        throw new IllegalArgumentException(
                "ARGModel.Node can only have two parents");
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
        return nodes.get(i);
    }

    public final NodeRef getInternalNode(int i) {
        return nodes.get(i + externalNodeCount);
    }

    public final NodeRef getNode(int i) {
        return nodes.get(i);
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


    public final int getReassortmentNodeCount() {
        int cnt = 0;
        for (Node node : nodes) {
            if (!node.bifurcation)
                cnt++;
        }
        return cnt;
    }

    // public final int getReassortmentNodeCount() { return nullCounter; }

    public void addNullCounter() {
        nullCounter++;
    }

    public void removeNullCounter() {
        nullCounter--;
    }

    /**
     * Returns the root node of this tree.
     */
    public final NodeRef getRoot() {
        return root;
    }

    public final NodeRef getRoot(int partition) {
        // TODO
        return null;
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    public final void setRoot(NodeRef newRoot) {

        if (!inEdit)
            throw new RuntimeException(
                    "Must be in edit transaction to call this method!");

        root = (Node) newRoot;

        // We shouldn't need this because the addChild will already have fired
        // appropriate events.
        // pushTreeChangedEvent();
    }

    public void swapHeightParameters(NodeRef n1, NodeRef n2) {
        Node node1 = (Node) n1;
        Node node2 = (Node) n2;
        double height1 = node1.getHeight();
        double height2 = node2.getHeight();
        Parameter trans = node1.heightParameter;
        node1.heightParameter = node2.heightParameter;
        node2.heightParameter = trans;
        node1.setHeight(height1);
        node2.setHeight(height2);

    }

    /**
     * Links <code>parent</code> with <code>child</code>.  If
     * <code>parent</code> is a bifurcation node,
     * the method calls <code>singleAddChild(parent,child)</code>,
     * otherwise, the method calls <code>doubleAddChild(parent,child)</code>.
     *
     * @throws RuntimeException if not in edit mode.
     * @see <code>singleAddChild(NodeRef parent, NodeRef child)</code>
     * @see <code>doubleAddChild(NodeRef parent, NodeRef child)</code>
     */
    public void addChild(NodeRef parent, NodeRef child) {
        checkEditMode();

        Node p = (Node) parent;
        Node c = (Node) child;

        if (p.bifurcation) {
            p.singleAddChild(c);
        } else {
            p.doubleAddChild(c);
        }
    }

    public void addChildWithSingleParent(NodeRef parent, NodeRef child) {
        checkEditMode();

        Node p = (Node) parent;
        Node c = (Node) child;

        if (p.bifurcation) {
            p.singleAddChildWithOneParent(c);
        } else {
            p.doubleAddChildWithOneParent(c);
        }

    }


    /**
     * Makes a link between <code>parent</code> and <code>child</code>.
     * By default, if <code>parent</code> has a null reference for both
     * it's children, <code>child</code> will become the left child of
     * <code>parent</code>, otherwise <code>child</code> will become
     * the right child of <code>parent</code>.  <br><br>If the right parent
     * of <code>child</code> is <code>null</code>, <code>parent</code>
     * will become the parent, the same thing will happen for the left parent
     * of <code>child</code>.
     *
     * @param parent the <code>NodeRef</code> that will become the parent of <code>child</code>
     * @param child  the <code>NodeRef</code> that will become the child of <code>parent</code>
     * @throws RuntimeException         if the you are not in edit transaction mode
     * @throws IllegalArgumentException if <code>parent</code> already has two children.
     */
    public void singleAddChild(NodeRef parent, NodeRef child) {

        if (!inEdit) {
            throw new RuntimeException(
                    "must be in edit transaction to call this method!");
        }

        Node p = (Node) parent;
        Node c = (Node) child;
        p.singleAddChild(c);
    }

    public void singleAddChildWithOneParent(NodeRef p, NodeRef c) {

        if (!inEdit)
            throw new RuntimeException(
                    "Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.singleAddChildWithOneParent(child);
    }

    public void doubleAddChild(NodeRef p, NodeRef c) {

        if (!inEdit)
            throw new RuntimeException(
                    "Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.doubleAddChild(child);
    }

    public void doubleAddChildWithOneParent(NodeRef p, NodeRef c) {

        if (!inEdit)
            throw new RuntimeException(
                    "Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.doubleAddChildWithOneParent(child);
    }

    public void addChildAsRecombinant(NodeRef p1, NodeRef p2, NodeRef c,
                                      Parameter partitioning) {
        // public void addChildAsRecombinant(NodeRef p1, NodeRef p2, NodeRef c,
        // BitSet bs1, BitSet bs2) {
        if (!inEdit)
            throw new RuntimeException(
                    "Must be in edit transaction to call this method!");
        Node parent1 = (Node) p1;
        Node parent2 = (Node) p2;
        Node child = (Node) c;
        // if (parent1.hasChild(child) || parent2.hasChild(child)) throw new
        // IllegalArgumentException("Child already exists in parent");
        // if (parent2.hasChild(child)) throw new
        // IllegalArgumentException("Child already exists in")
        parent1.addChildRecombinant(child, partitioning);
        parent2.addChildRecombinant(child, partitioning);
        // if( parent2.getChildCount() == 1 )
        // parent2.addChildNoParentConnection(node);
    }

    /**
     * Removes the link between <code>parent</code> and
     * <code>child</code>.  If <code>parent</code> is a bifurcation node,
     * the method calls <code>singleRemoveChild(parent, child)</code>,
     * otherwise, the method calls <code>doubleRemoveChild(parent,child)</code>.
     *
     * @throws RuntimeException if not in edit mode.
     * @see <code>singleRemoveChild(NodeRef parent, NodeRef child)</code>
     * @see <code>doubleRemoveChild(NodeRef parent, NodeRef child)</code>
     */
    public void removeChild(NodeRef parent, NodeRef child) {
        checkEditMode();

        Node p = (Node) parent;
        Node c = (Node) child;

        p.doubleRemoveChild(c);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {

    }

    /**
     * Removes the link between the parent and the child. This method should be
     * called when the child's parents are both the same. After the method is
     * called the parent will have two null references for children, and the
     * child will have two null references for parents.
     *
     * @param parent The parent NodeRef
     * @param child  The child NodeRef
     *               //	 * @see doubleAddChild()
     *               //	 * @see singleAddChild()
     *               //	 * @see removeChild()
     *               //	 * @see singleRemoveChild()
     */
    public void doubleRemoveChild(NodeRef parent, NodeRef child) {
        checkEditMode();

        Node p = (Node) parent;
        Node c = (Node) child;

        p.doubleRemoveChild(c);
    }

    public void singleRemoveChild(NodeRef p, NodeRef c) {
        checkEditMode();

        Node parent = (Node) p;
        Node child = (Node) c;

        parent.singleRemoveChild(child);
    }

    protected Node oldRoot;


    public boolean beginTreeEdit() {
        if (inEdit)
            throw new RuntimeException("Alreading in edit transaction mode!");

        oldRoot = root;

        inEdit = true;

        return false;
    }

    public void endTreeEdit() {
        if (!inEdit)
            throw new RuntimeException("Not in edit transaction mode!");

        inEdit = false;

        if (root != oldRoot) {
            swapParameterObjects(oldRoot, root);
        }

        // ystem.err.println("There are "+treeChangedEvents.size()+" events
        // waiting");
        // System.exit(-1);
        for(TreeChangedEvent treeChangedEvent : treeChangedEvents) {
            listenerHelper.fireModelChanged(this, treeChangedEvent);
        }
        treeChangedEvents.clear();
    }

    public void checkTreeIsValid() throws MutableTree.InvalidTreeException {
        for (Node node : nodes) {
            if (!node.heightParameter.isWithinBounds()) {
                throw new InvalidTreeException("height parameter out of bounds");
            }
        }
    }

    private void endTreeEditFast() {
        inEdit = false;
    }

    public void setNodeHeight(NodeRef n, double height) {
        ((Node) n).setHeight(height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        if (!hasRates)
            throw new IllegalArgumentException(
                    "Rate parameters have not been created");
        ((Node) n).setRate(rate);

    }

    public void setNodeTrait(NodeRef n, double value) {
        if (!hasTraits)
            throw new IllegalArgumentException(
                    "Trait parameters have not been created");
        ((Node) n).setTrait(value);
    }

    public void setNodeNumber(NodeRef node, int n) {
        node.setNumber(n);
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException(
                "ARGModel cannot have branch lengths set");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException(
                "ARGModel does not use NodeAttributes");
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    /**
     * Store current state
     */
    protected void storeState() {
        /*
           * System.err.println("Storing state"); this.checkBranchSanity();
           * System.err.println("sane before operation");
           */
        copyNodeStructure(storedNodes);
        // storedRootNumber = storedNodes.indexOf(root.getNumber();
        storedRootNumber = nodes.indexOf(root);
        storedNodeCount = nodeCount;
        storedInternalNodeCount = internalNodeCount;
        addedParameters = null;
        addedPartitioningParameter = null;
        addedNodes = null;
        removedParameters = null;
        removedPartitioningParameter = null;
        removedNodes = null;
        storedNullCounter = nullCounter;
        // System.err.println("Stored: "+Tree.Utils.uniqueNewick(this,
        // getRoot()));
        // System.err.println("Stored : "+this.toString());
    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        ArrayList<Node> tmp = storedNodes;
        storedNodes = nodes;
        nodes = tmp;
        root = nodes.get(storedRootNumber);

        nodeCount = storedNodeCount;
        internalNodeCount = storedInternalNodeCount;

        //This part occurs true when we add a new arg event
        //onto the root.

        if (addedParameters != null) {
            if (addedParameters.length == 5) {
                removeVariable(addedParameters[0]);
                removeVariable(addedParameters[1]);
                removeVariable(addedParameters[2]);
                removeVariable(addedParameters[3]);

                storedInternalNodeHeights.removeParameter(addedParameters[0]);
                storedInternalNodeHeights.removeParameter(addedParameters[4]);

                storedInternalAndRootNodeHeights.removeParameter(addedParameters[0]);
                storedInternalAndRootNodeHeights.removeParameter(addedParameters[1]);

                storedNodeRates.removeParameter(addedParameters[2]);
            } else {
                storedInternalNodeHeights.removeParameter(addedParameters[0]);
                storedInternalNodeHeights.removeParameter(addedParameters[1]);

                removeVariable(addedParameters[0]);
                removeVariable(addedParameters[1]);

                storedInternalAndRootNodeHeights
                        .removeParameter(addedParameters[0]);
                storedInternalAndRootNodeHeights
                        .removeParameter(addedParameters[1]);

                removeVariable(addedParameters[2]);
                removeVariable(addedParameters[3]);
                storedNodeRates.removeParameter(addedParameters[2]);
//                storedNodeRates.removeVariable(addedParameters[3]);
            }
        }
        if (addedPartitioningParameter != null) {
            partitioningParameters.removeParameter(addedPartitioningParameter);
            removeVariable(addedPartitioningParameter);
        }
        if (removedParameters != null) {
            storedInternalNodeHeights.addParameter(removedParameters[0]);
            storedInternalNodeHeights.addParameter(removedParameters[1]);
            addVariable(removedParameters[0]);
            addVariable(removedParameters[1]);

            storedInternalAndRootNodeHeights.addParameter(removedParameters[0]);
            storedInternalAndRootNodeHeights.addParameter(removedParameters[1]);

            addVariable(removedParameters[2]);
            addVariable(removedParameters[3]);
            storedNodeRates.addParameter(removedParameters[2]);
//            storedNodeRates.addVariable(removedParameters[3]);

        }

        if (removedPartitioningParameter != null) {
            partitioningParameters.addParameter(removedPartitioningParameter);
            addVariable(removedPartitioningParameter);
        }


        nullCounter = storedNullCounter;
    }

    /**
     * accept the stored state
     */
    protected void acceptState() {
        // System.err.println("Accepted ARG\n"+this.toGraphString());
    } // nothing to do

    /**
     * Adopt the state of the model component from source.
     */
    protected void adoptState(Model source) {
    }

    /*
      * public void addNewHeightParameters(Parameter newbie1, Parameter newbie2,
      * CompoundParameter internalNodeParameters,
      * CompoundParameter internalAndRootNodeParameters) {
      * addVariable(newbie1); addVariable(newbie2); addedParameters = new
      * Parameter[2]; addedParameters[0] = newbie1; addedParameters[1] = newbie2;
      *
      * storedInternalNodeHeights = internalNodeParameters;
      * storedInternalNodeHeights.addVariable(newbie1);
      * storedInternalNodeHeights.addVariable(newbie2);
      *
      *//*
		 * storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
		 * storedInternalAndRootNodeHeights.addVariable(newbie1);
		 * storedInternalAndRootNodeHeights.addVariable(newbie2);
		 *//*
																 * }
																 */

    public void expandARG(Node newbie1, Node newbie2,
                          CompoundParameter internalNodeParameters,
                          CompoundParameter internalAndRootNodeParameters,
                          CompoundParameter nodeRates) {

        addVariable(newbie1.heightParameter);
        addVariable(newbie2.heightParameter);
        addVariable(newbie2.partitioning);

        addVariable(newbie1.rateParameter);
        addVariable(newbie2.rateParameter);

        addedParameters = new Parameter[4];
        addedParameters[0] = newbie1.heightParameter;
        addedParameters[1] = newbie2.heightParameter;
        addedParameters[2] = newbie1.rateParameter;
        addedParameters[3] = newbie2.rateParameter;
        addedPartitioningParameter = newbie2.partitioning;

        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalNodeHeights.addParameter(newbie2.heightParameter);

        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalAndRootNodeHeights.addParameter(newbie2.heightParameter);

        storedNodeRates = nodeRates;
        storedNodeRates.addParameter(newbie1.rateParameter);
//         storedNodeRates.addVariable(newbie2.rateParameter);


        partitioningParameters.addParameter(newbie2.partitioning);
        nodes.add(newbie1);
        nodes.add(newbie2);
        internalNodeCount += 2;

//         pushTreeSizeIncreasedEvent();

    }

    public void expandARGWithRecombinant(Node newbie1, Node newbie2,
                                         CompoundParameter internalNodeParameters,
                                         CompoundParameter internalAndRootNodeParameters,
                                         CompoundParameter nodeRates) {
        // System.err.println("attempting to expand");

        addVariable(newbie1.heightParameter);
        addVariable(newbie2.heightParameter);
        addVariable(newbie2.partitioning);

        // System.err.println("expand 0");

        addVariable(newbie1.rateParameter);
        addVariable(newbie2.rateParameter);

        // System.err.println("expand 1");

        addedParameters = new Parameter[4];
        addedParameters[0] = newbie1.heightParameter;
        addedParameters[1] = newbie2.heightParameter;
        addedParameters[2] = newbie1.rateParameter;
        addedParameters[3] = newbie2.rateParameter;
        addedPartitioningParameter = newbie2.partitioning;

        // System.err.println("expand 2");

        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalNodeHeights.addParameter(newbie2.heightParameter);

        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalAndRootNodeHeights.addParameter(newbie2.heightParameter);

        storedNodeRates = nodeRates;
        storedNodeRates.addParameter(newbie1.rateParameter);
        storedNodeRates.addParameter(newbie2.rateParameter);

        // System.err.println("expand 3");

        partitioningParameters.addParameter(newbie2.partitioning);
        nodes.add(newbie1);
        nodes.add(newbie2);
        internalNodeCount += 2;
        // sanityNodeCheck(internalNodeParameters);
//        pushTreeSizeIncreasedEvent();

        // System.err.println("done expand");

    }

    public void sanityNodeCheck(CompoundParameter inodes) {
        int len = inodes.getParameterCount();
        for (int i = 0; i < len; i++) {
            Parameter p = inodes.getParameter(i);
            for (int j = 0; j < internalNodeCount; j++) {
                Node node = (Node) getInternalNode(j);
                if (node.heightParameter == p) {
                    if (isRoot(node)) {
                        System.err
                                .println("Root height found in internal nodes");
                        System.exit(-1);
                    }
                }

            }
        }
    }

    public void contractARGWithRecombinantNewRoot(Node oldie,
                                                  Node oldRoot, Node newRoot,
                                                  CompoundParameter internalNodeParameters,
                                                  CompoundParameter internalAndRootNodeParameters,
                                                  CompoundParameter nodeRates) {
        removeVariable(oldie.heightParameter);
        removeVariable(oldRoot.heightParameter);
        removeVariable(oldie.partitioning);

        removeVariable(oldie.rateParameter);
        removeVariable(oldie.rateParameter);

        removedParameters = new Parameter[4];
        removedParameters[0] = oldie.heightParameter;
        removedParameters[1] = oldRoot.heightParameter;
        removedParameters[2] = oldie.rateParameter;
        removedParameters[3] = oldRoot.rateParameter;

        partitioningParameters.removeParameter(oldie.partitioning);
        removedPartitioningParameter = oldie.partitioning;
        storedInternalNodeHeights = internalNodeParameters;

        storedInternalNodeHeights.removeParameter(oldie.heightParameter);

        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights
                .removeParameter(oldie.heightParameter);
        storedInternalAndRootNodeHeights
                .removeParameter(oldRoot.heightParameter);

        storedNodeRates = nodeRates;
        storedNodeRates.removeParameter(oldie.rateParameter);
        storedNodeRates.removeParameter(oldRoot.rateParameter);

        nodes.remove(oldie);
        nodes.remove(oldRoot);
        internalNodeCount -= 2;

        this.setRoot(newRoot);
//        pushTreeSizeDecreasedEvent();


    }


    public void contractARG(Node oldie1, Node oldie2,
                            CompoundParameter internalNodeParameters,
                            CompoundParameter internalAndRootNodeParameters,
                            CompoundParameter nodeRates) {

        removeVariable(oldie1.heightParameter);
        removeVariable(oldie2.heightParameter);
        removeVariable(oldie2.partitioning);


        removeVariable(oldie1.rateParameter);
        removeVariable(oldie2.rateParameter);


        removedParameters = new Parameter[4];
        removedParameters[0] = oldie1.heightParameter;
        removedParameters[1] = oldie2.heightParameter;
        removedParameters[2] = oldie1.rateParameter;
        removedParameters[3] = oldie2.rateParameter;


        partitioningParameters.removeParameter(oldie2.partitioning);
        removedPartitioningParameter = oldie2.partitioning;
        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.removeParameter(oldie1.heightParameter);
        storedInternalNodeHeights.removeParameter(oldie2.heightParameter);


        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights.removeParameter(oldie1.heightParameter);
        storedInternalAndRootNodeHeights.removeParameter(oldie2.heightParameter);


        storedNodeRates = nodeRates;
        storedNodeRates.removeParameter(oldie1.rateParameter);
//    	storedNodeRates.removeVariable(oldie2.rateParameter);

        nodes.remove(oldie1);
        nodes.remove(oldie2);

        internalNodeCount -= 2;
//    	pushTreeSizeDecreasedEvent();
    }


    /**
     * Cleans up the arg model after a deletion event.
     *
     * @param oldie1
     * @param oldie2
     * @param internalNodeParameters
     * @param internalAndRootNodeParameters
     * @param nodeRates
     */
    public void contractARGWithRecombinant(Node oldie1, Node oldie2,
                                           CompoundParameter internalNodeParameters,
                                           CompoundParameter internalAndRootNodeParameters,
                                           CompoundParameter nodeRates) {

        removeVariable(oldie1.heightParameter);
        removeVariable(oldie2.heightParameter);
        removeVariable(oldie2.partitioning);


        removeVariable(oldie1.rateParameter);
        removeVariable(oldie2.rateParameter);


        removedParameters = new Parameter[4];
        removedParameters[0] = oldie1.heightParameter;
        removedParameters[1] = oldie2.heightParameter;
        removedParameters[2] = oldie1.rateParameter;
        removedParameters[3] = oldie2.rateParameter;


        partitioningParameters.removeParameter(oldie2.partitioning);
        removedPartitioningParameter = oldie2.partitioning;
        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.removeParameter(oldie1.heightParameter);
        storedInternalNodeHeights.removeParameter(oldie2.heightParameter);


        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights
                .removeParameter(oldie1.heightParameter);
        storedInternalAndRootNodeHeights
                .removeParameter(oldie2.heightParameter);


        storedNodeRates = nodeRates;
        storedNodeRates.removeParameter(oldie1.rateParameter);
        storedNodeRates.removeParameter(oldie2.rateParameter);

        nodes.remove(oldie1);
        nodes.remove(oldie2);

        internalNodeCount -= 2;
//        pushTreeSizeDecreasedEvent();
    }

    public boolean argStoreCheck() {
        if (storedInternalNodeHeights.getDimension() == this.internalNodeCount)
            return true;
        return false;
    }

    /**
     * Copies the node connections from this ARGModel's nodes array to the
     * destination array. Basically it connects up the nodes in destination in
     * the same way as this ARGModel is set up. This method is package private.
     */
    void copyNodeStructure(ArrayList<Node> destination) {

        // if ( nodes.length != destination.length ) {
        // throw new IllegalArgumentException("Node arrays are of different
        // lengths");
        // }
        while (destination.size() < nodes.size())
            destination.add(new Node());
        while (destination.size() > nodes.size())
            destination.remove(0);
        int n = nodes.size();
        // System.err.println("node.size = "+n);
        for (int i = 0; i < n; i++) {
            // for( Node node0 : nodes ) {
            Node node0 = nodes.get(i);
            Node node1 = destination.get(i);

            // the parameter values are automatically stored and restored
            // just need to keep the links
            node1.heightParameter = node0.heightParameter;
            node1.rateParameter = node0.rateParameter;
            node1.traitParameter = node0.traitParameter;
            node1.partitioning = node0.partitioning;

            node1.taxon = node0.taxon;
            node1.bifurcation = node0.bifurcation;
            node1.number = node0.number;
            node1.myHashCode = node0.myHashCode;
            // node1.partitionSet = (BitSet)node0.partitionSet.clone();
            // if (node0.leftPartition != null) {
            // node1.leftPartition = (BitSet) node0.leftPartition.clone();
            // } else {
            // node1.leftPartition = null;
            // }
            // if (node0.rightPartition != null) {
            // node1.rightPartition = (BitSet) node0.rightPartition.clone();
            // } else {
            // node1.rightPartition = null;
            // }
            //

            if (node0.leftParent != null) {
                node1.leftParent = // storedNodes.get(node0.leftParent.getNumber());
                        storedNodes.get(nodes.indexOf(node0.leftParent));
            } else {
                node1.leftParent = null;
            }

            if (node0.rightParent != null) {
                node1.rightParent = // storedNodes.get(node0.rightParent.getNumber());
                        storedNodes.get(nodes.indexOf(node0.rightParent));
            } else {
                node1.rightParent = null;
            }

            if (node0.leftChild != null) {
                node1.leftChild = // storedNodes.get(node0.leftChild.getNumber());
                        storedNodes.get(nodes.indexOf(node0.leftChild));
            } else {
                node1.leftChild = null;
            }

            if (node0.rightChild != null) {
                node1.rightChild = // storedNodes.get(node0.rightChild.getNumber());
                        storedNodes.get(nodes.indexOf(node0.rightChild));
            } else {
                node1.rightChild = null;
            }
        }
    }

    public void setPartitionRecursively(NodeRef nr, int partition) {
        Node node = (Node) nr;
        node.setPartitionRecursively(partition);
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
        if (i == 0)
            return root.heightParameter;
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
            if (getTaxonId(i).equals(id))
                return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon)
                return i;
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
                index++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the taxon of the
     *         given external node. If the node doesn't have a taxon then the
     *         nodes own attribute is returned.
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
        throw new IllegalArgumentException("Cannot add taxon to a ARGModel");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a ARGModel");
    }

    public void setTaxonId(int taxonIndex, String id) {
        throw new IllegalArgumentException("Cannot set taxon id in a ARGModel");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new IllegalArgumentException(
                "Cannot set taxon attribute in a ARGModel");
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
    public final String getNewick(int partition) {
        return TreeUtils.newick(new ARGTree(this, partition));
        // return Tree.Utils.newick(this);
    }

    /**
     * Checks whether <code>ARGMode</code> is in edit mode.
     *
     * @throws RuntimeException if the <code>ARGModel</code> is not in edit mode.
     */
    private void checkEditMode() throws RuntimeException {
        if (!inEdit)
            throw new RuntimeException("Not in edit transaction mode!");
    }

    public void checkBranchSanity() {
        boolean plotted = false;
        for (Node node : nodes) {
            if (!node.isRoot()) {
                double length1 = 0;
                double length2 = 0;
                if (node.leftParent != null)
                    length1 = getNodeHeight(node.leftParent)
                            - getNodeHeight(node);
                if (node.rightParent != null)
                    length2 = getNodeHeight(node.rightParent)
                            - getNodeHeight(node);
                if (String.valueOf(length1).equals("NaN")
                        || String.valueOf(length2).equals("NaN")) {
                    if (!plotted) {
                        System.err.println(toGraphString());
                        plotted = true;
                    }
                    System.err.println("Caught the NaN: node=" + node.number
                            + " (" + node.getHeight() + ") lp="
                            + node.leftParent.number + " ("
                            + node.leftParent.getHeight() + ") rp="
                            + node.rightParent.number + " ("
                            + node.rightParent.getHeight() + ")");
                    System.exit(-1);
                }
            }
        }

    }

    /**
     * @return a string containing an extended newick representation of the tree
     */
    public String toString() {
        // StringBuffer sb = new StringBuffer();
        // for (int i = 0; i < maxNumberOfPartitions; i++) {
        // sb.append(i + ": ");
        // sb.append(getNewick(i));
        // }
        // return new String(sb);
        return toExtendedNewick();
    }

    public static final String nullEdge = " -";

    public void appendGraphStringOld(StringBuffer sb) {
        int cnt = 0;

        for (Node node : nodes)
            node.number = cnt++;

        cnt = 0;
        for (Node node : nodes) {
            sb.append(cnt == 0 ? "[" : ",[");
            cnt++;
            sb.append(node.number + ":");

            if (node.leftParent == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.leftParent.number);

            if (node.rightParent == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.rightParent.number);

            if (node.leftChild == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.leftChild.number);

            if (node.rightChild == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.rightChild.number);
            // sb.append(" " + node.bifurcation);
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
            // if (node.leftPartition != null)
            // sb.append(" l");
            // if (node.rightPartition != null)
            // sb.append(" r");
            sb.append("]");
            //
            // );
        }
        // sb.append("Root = " + ((Node) getRoot()).number + "\n");
        // sb.append("\n");
    }

    public boolean validRoot() {
        // todo -- there must be a way to some graph properties to do this
        // check.
        boolean valid = true;
        for (int i = 0; valid && i < maxNumberOfPartitions; i++) {
            ARGTree argTree = new ARGTree(this, i);
            if (argTree.wasRootTrimmed())
                valid = false;
        }
        return valid;
    }

    public String toGraphString() {
        int cnt = 1;
        for (Node node : nodes) {
            node.number = cnt;
            cnt++;
        }
        StringBuffer sb = new StringBuffer();
        sb.append("Total length: " + nodes.size() + "\n");
        for (Node node : nodes) {
            sb.append(node.number + ":");

            if (node.leftParent == null)
                sb.append(" 0");
            else
                sb.append(" " + node.leftParent.number);

            if (node.rightParent == null)
                sb.append(" 0");
            else
                sb.append(" " + node.rightParent.number);

            if (node.leftChild == null)
                sb.append(" 0");
            else
                sb.append(" " + node.leftChild.number);

            if (node.rightChild == null)
                sb.append(" 0");
            else
                sb.append(" " + node.rightChild.number);
            // sb.append(" " + node.bifurcation);
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
            if (node.partitioning != null)
                sb.append(" p");
            /*
                * if (node.leftPartition != null) sb.append(" l"); if
                * (node.rightPartition != null) sb.append(" r");
                */
            sb.append("\t" + getNodeHeight(node));
            sb.append("\n");
        }
        sb.append("Root = " + ((Node) getRoot()).number + "\n");
        return new String(sb);
    }

    public ARGModel fromGraphStringCompressed(String source) {
        StringTokenizer st1 = new StringTokenizer(source, ":");
        int numberNodes = Integer.parseInt(st1.nextToken());
        int numberPartitions = Integer.parseInt(st1.nextToken());
        int rootNumber = Integer.parseInt(st1.nextToken());
        int external = 0;

        ArrayList<Node> nodes = new ArrayList<Node>();
        for (int i = 0; i < numberNodes; i++) {
            Node node = new Node();
            node.number = i;
            nodes.add(node);
        }

        for (int i = 0; i < numberNodes; i++) {
            String nodeString = st1.nextToken();
            StringTokenizer st2 = new StringTokenizer(nodeString);
            Node node = nodes.get(i);
            int lP = Integer.parseInt(st2.nextToken());
            int rP = Integer.parseInt(st2.nextToken());
            int lC = Integer.parseInt(st2.nextToken());
            int rC = Integer.parseInt(st2.nextToken());
            if (lP != -1)
                node.leftParent = nodes.get(lP);
            if (rP != -1)
                node.rightParent = nodes.get(rP);
            if (lC != -1)
                node.leftChild = nodes.get(lC);
            if (rC != -1)
                node.rightChild = nodes.get(rC);
            double height = Double.parseDouble(st2.nextToken());
            node.heightParameter = new Parameter.Default(height);
            String name = st2.nextToken();
            if (name.compareTo("NA") != 0) {
                node.taxon = new Taxon(name); // todo saveList
                external++;
            }
            String p = st2.nextToken();
            if (p.compareTo("NA") == 0)
                node.bifurcation = true;
            else {
                node.bifurcation = false;
                node.partitioning = new Parameter.Default(numberPartitions, 0.0);
                node.partitioning.setParameterValueQuietly(Integer.parseInt(p),
                        1.0);
                while (st2.hasMoreTokens())
                    node.partitioning.setParameterValueQuietly(Integer
                            .parseInt(st2.nextToken()), 1.0);
            }
        }

        return new ARGModel(nodes, nodes.get(rootNumber), numberPartitions,
                external);

        // return null;
    }


    /**
     * Gives a summary of the ARG model.  Should only
     * be used for debugging purposes because it's really slow.
     *
     * @return a summary of the ARG model.
     */
    public String toARGSummary() {
        NumberFormatter format = new NumberFormatter(4);

        String space = "   ";

        String a = "----------------------\n" +
                "ARG Summary \n---------------------- \n";

        a += "Number of nodes: " + nodes.size() + "\n";
        a += "Number of partitions: " + maxNumberOfPartitions + "\n";
        a += "Number of Reassorments: " + this.getReassortmentNodeCount() + "\n";
        a += "Root number: " + getRoot().getNumber() + "\n";
        a += "Node Summary"
                + "\n----------------------------------------\n" +
                "ID  LP   RP   LC   RC   Height" + space + "TX  \n" +
                "----------------------------------------\n";

        for (Node node : nodes) {
            a += node.getNumber() + space;

            if (node.leftParent == null) a += "-1" + space;
            else a += " " + node.leftParent.number + space;
            if (node.rightParent == null) a += "-1" + space;
            else a += " " + node.rightParent.number + space;
            if (node.leftChild == null) a += "-1" + space;
            else a += " " + node.leftChild.number + space;
            if (node.rightChild == null) a += "-1" + space;
            else a += " " + node.rightChild.number + space;

            a += format.formatDecimal(getNodeHeight(node), 4) + space;

            if (node.partitioning != null) {
                for (int i = 0, n = getNumberOfPartitions(); i < n; i++) {
                    a += node.partitioning.getParameterValue(i) + space;
                }
            }

            if (node.taxon == null) {
                a += "internal" + space;
            } else {
                a += node.taxon + space;
            }
            a += "\n";
        }
        a += "\nInduced Trees" +
                "\n----------------------------------------\n";
//		for(int i = 0; i < maxNumberOfPartitions; i++){
//			ARGTree tree = new ARGTree(this,i);
//			a += "Partition " + i + "\n  " + tree.toString() + "\n";
//		}

        return a;
    }

    public String toGraphStringCompressed(boolean recurse) {
        int cnt = 0;
        for (Node node : nodes) {
            node.number = cnt;
            cnt++;
        }
        StringBuffer sb = new StringBuffer();
        // sb.append("Total length: " + nodes.size() + "\n");
        sb.append(nodes.size());
        sb.append(":");
        sb.append(maxNumberOfPartitions);
        sb.append(":");
        sb.append(getRoot().getNumber());
        for (Node node : nodes) {
            sb.append(":");

            if (node.leftParent == null)
                sb.append(" -1");
            else
                sb.append(" " + node.leftParent.number);

            if (node.rightParent == null)
                sb.append(" -1");
            else
                sb.append(" " + node.rightParent.number);

            if (node.leftChild == null)
                sb.append(" -1");
            else
                sb.append(" " + node.leftChild.number);

            if (node.rightChild == null)
                sb.append(" -1");
            else
                sb.append(" " + node.rightChild.number);
            sb.append(" " + getNodeHeight(node));
            // sb.append(" " + node.bifurcation);
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
            else
                sb.append(" NA");
            if (node.partitioning != null) {
                // sb.append(" p");
                double[] bits = node.partitioning.getParameterValues();
                for (int i = 0; i < bits.length; i++) {
                    if (bits[i] == 1) {
                        sb.append(" " + i);
                    }
                }

            } else {
                sb.append(" NA");
            }

            /*
                * if (node.leftPartition != null) sb.append(" l"); if
                * (node.rightPartition != null) sb.append(" r");
                */
            // sb.append("\t" + getNodeHeight(node));
            // sb.append("\n");
        }
        // sb.append("Root = " + ((Node) getRoot()).number + "\n");

        String rtn = new String(sb);

        /*
           * if ( recurse ) { ARGModel test = fromGraphStringCompressed(rtn);
           * System.err.println("OLD 0: "+getNewick(0)); System.err.println("NEW
           * 0: "+test.getNewick(0)); System.err.println("OLD:
           * "+toGraphStringCompressed(false)); System.err.println("NEW:
           * "+test.toGraphStringCompressed(false)); }
           */

        return new String(sb);
    }

    public Tree getCopy() {
        throw new UnsupportedOperationException(
                "please don't call this function");
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
    protected Node getNodeOfParameter(Parameter parameter) {

        if (parameter == null)
            throw new IllegalArgumentException("Parameter is null!");

        // for (int i =0; i < nodes.length; i++) {
        for (Node node : nodes) {
            if (node.heightParameter == parameter) {
                return node;
            }
            if (hasRates && node.rateParameter == parameter) {
                return node;
            }
            if (hasTraits && node.traitParameter == parameter) {
                return node;
            }
        }
        throw new RuntimeException("Parameter not found in any nodes:"
                + parameter.getId());
    }

    /**
     * Get the root height parameter. Is private because it can only be called
     * by the XMLParser
     */
    public Parameter getRootHeightParameter() {

        return root.heightParameter;
    }

    /**
     * @return the relevant node height parameter. Is private because it can
     *         only be called by the XMLParser
     */
    public Parameter createNodeHeightsParameter(boolean rootNode,
                                                boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException(
                    "At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter("nodeHeights");
        // System.err.println("Constructed nodeHeights");
        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.heightParameter);
            }
        }

        if (leafNodes) {
            for (int i = 0; i < externalNodeCount; i++) {
                parameter.addParameter(nodes.get(i).heightParameter);
            }
        }

        return parameter;
    }

    public Parameter getLeafHeightParameter(NodeRef nr) {
        Node node = (Node) nr;
        if (!isExternal(node)) {
            throw new RuntimeException(
                    "only root and leaves can be used with setNodeHeightParameter");
        }

        return nodes.get(nodes.indexOf(node)).heightParameter;
    }

    /**
     * @return the relevant node rate parameter. Is private because it can only
     *         be called by the XMLParser
     */
    public Parameter createNodeRatesParameter(boolean rootNode,
                                              boolean internalNodes, boolean leafNodes, int numberPartitions) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException(
                    "At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter(TreeModelParser.NODE_RATES);

        hasRates = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            node.createRateParameter(numberPartitions);
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.rateParameter);
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            Node node = nodes.get(i);
            node.createRateParameter(numberPartitions);
            if (leafNodes) {
                parameter.addParameter(node.rateParameter);
            }
        }

        return parameter;
    }

    /**
     * Create a node traits parameter. Is private because it can only be called
     * by the XMLParser
     */
    public Parameter createNodeTraitsParameter(boolean rootNode,
                                               boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException(
                    "At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter(TreeModelParser.NODE_TRAITS);

        hasTraits = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            node.createTraitParameter();
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.traitParameter);
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            Node node = nodes.get(i);
            node.createTraitParameter();
            if (leafNodes) {
                parameter.addParameter(node.traitParameter);
            }
        }

        return parameter;
    }

    /**
     * This method swaps the parameter objects of the two nodes but maintains
     * the values in each node. This method is used to ensure that root node of
     * the tree always has the same parameter object.
     */
    private void swapParameterObjects(Node n1, Node n2) {

        double height1 = n1.getHeight();
        double height2 = n2.getHeight();

        double rate1 = 1.0, rate2 = 1.0;
        double trait1 = 0.0, trait2 = 0.0;

        if (hasRates) {
            System.exit(-1);

            rate1 = n1.getRate(0);
            rate2 = n2.getRate(0);
        }

        if (hasTraits) {
            trait1 = n1.getTrait();
            trait2 = n2.getTrait();
        }

        Parameter temp = n1.heightParameter;
        n1.heightParameter = n2.heightParameter;
        n2.heightParameter = temp;

        if (hasRates) {
            temp = n1.rateParameter;
            n1.rateParameter = n2.rateParameter;
            n2.rateParameter = temp;
        }

        if (hasTraits) {
            temp = n1.traitParameter;
            n1.traitParameter = n2.traitParameter;
            n2.traitParameter = temp;
        }

        n1.heightParameter.setParameterValueQuietly(0, height1);
        n2.heightParameter.setParameterValueQuietly(0, height2);

        if (hasRates) {
            n1.rateParameter.setParameterValueQuietly(0, rate1);
            n2.rateParameter.setParameterValueQuietly(0, rate2);
        }

        if (hasTraits) {
            n1.traitParameter.setParameterValueQuietly(0, trait1);
            n2.traitParameter.setParameterValueQuietly(0, trait2);
        }

    }

    // **************************************************************
    // Private inner classes
    // **************************************************************
    public class Node implements NodeRef {

        public boolean[] ancestralMaterial;
        public boolean fullAncestralMaterial;
        public boolean hasSomeAncestralMaterial;

        public boolean hasReassortmentAncestor() {
            Node a = this;

            while (a != null) {
                if (!a.bifurcation) {
                    return true;
                } else {
                    a = a.leftParent;
                }
            }


            return false;
        }

        public void setAncestralMaterial(boolean[] childAncestralMaterial) {
            if (fullAncestralMaterial) {
                return;
            }

            fullAncestralMaterial = true;

            for (int i = 0; i < ancestralMaterial.length; i++) {
                ancestralMaterial[i] = ancestralMaterial[i] || childAncestralMaterial[i];

                fullAncestralMaterial = fullAncestralMaterial && ancestralMaterial[i];
                hasSomeAncestralMaterial = hasSomeAncestralMaterial || ancestralMaterial[i];
            }

            if (bifurcation) {
                if (this.leftParent != null)
                    this.leftParent.setAncestralMaterial(ancestralMaterial);
            } else {
                boolean[] leftAncestralMaterial = new boolean[ancestralMaterial.length];
                boolean[] rightAncestralMaterial = new boolean[ancestralMaterial.length];

                System.arraycopy(ancestralMaterial, 0, leftAncestralMaterial, 0, leftAncestralMaterial.length);
                System.arraycopy(ancestralMaterial, 0, rightAncestralMaterial, 0, rightAncestralMaterial.length);

                for (int i = 0; i < ancestralMaterial.length; i++) {
                    if (partitioning.getParameterValue(i) == 0.0) {
                        rightAncestralMaterial[i] = false;
                    } else {
                        leftAncestralMaterial[i] = false;
                    }
                }

                this.leftParent.setAncestralMaterial(leftAncestralMaterial);
                this.rightParent.setAncestralMaterial(rightAncestralMaterial);
            }
        }


        public int myHashCode = 0;

        public int hashCode() {
            if (myHashCode == 0) {
                myHashCode = super.hashCode();
            }
            return myHashCode;
        }

        public boolean equals(Object o) {
            return hashCode() == o.hashCode();
        }


        public NodeRef mirrorNode;

        public Node leftParent, rightParent;

        public Node leftChild, rightChild;

        public int number;

        public Parameter heightParameter;

        //First half of the rate parameter represent the rates
        //Second half represents 0-1 indicators

        public Parameter rateParameter = null;

        public Parameter traitParameter = null;

        public Taxon taxon = null;

        // public BitSet leftPartition = null;
        // public BitSet rightPartition = null;
        // public Node dupSister = null;
        // public Node linkSister = null;
        // public Node dupParent = null;

        // public Node leftParent;
        // public Node rightParent;

        // public int leftPartition;
        // public int rightPartition;

        public Parameter partitioning;

        public boolean bifurcation = true;

        public int countReassortmentChild(Tree tree) {
            // int cnt = 0;
            if (isExternal())
                return 0;
            // if( leftChild == null ) {
            // System.err.println("left is null");
            // System.err.println("is reassort = "+isReassortment());
            // System.err.println("right child = "+Tree.Utils.uniqueNewick(tree,
            // rightChild));
            // }
            // if( rightChild == null ) {
            // System.err.println("right is null");
            // System.err.println("is reassort = "+isReassortment());
            // System.err.println("left child = "+Tree.Utils.uniqueNewick(tree,
            // leftChild));
            // }
            if (isReassortment()) {
                return 1 + leftChild.countReassortmentChild(tree);
            } else if (isBifurcationDoublyLinked()) {
                return leftChild.countReassortmentChild(tree);
            } else {
                return leftChild.countReassortmentChild(tree)
                        + rightChild.countReassortmentChild(tree);
            }
        }

        public Node() {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            heightParameter = null;
            number = 0;
            taxon = null;
            // leftPartition = null;
            // rightPartition = null;
            partitioning = null;
        }

        /**
         * Constructor to build an ARG into a bifurcating tree
         *
         * @param node
         */
        public Node(Node node) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            heightParameter = node.heightParameter;
            taxon = node.taxon;
            number = node.number;
            // if( node.isReassortment() ) {
            // Node parent = leftParent;
            // parent.removeChild(this);
            // return new Node(leftChild);
            // } else {
            if (node.leftChild != null) {
                if (node.leftChild.isReassortment())
                    singleAddChild(new Node(node.leftChild.leftChild));
                else
                    singleAddChild(new Node(node.leftChild));
            }
            if (node.rightChild != null) {
                if (node.rightChild.isReassortment())
                    singleAddChild(new Node(node.rightChild.leftChild));
                else
                    singleAddChild(new Node(node.rightChild));
            }
            // }
        }

        /**
         * constructor used to clone a node and all children with no
         * reassortments
         */
        public Node(Tree tree, NodeRef node) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            // leftPartition = new BitSet();
            // leftPartition.set(0);
            // leftPartition = null;
            // rightPartition = new BitSet();
            // rightPartition.set(0);
            // rightPartition = null;
            partitioning = null;
            heightParameter = new Parameter.Default(tree.getNodeHeight(node));
            addVariable(heightParameter);

            number = node.getNumber();
            taxon = tree.getNodeTaxon(node);

            for (int i = 0; i < tree.getChildCount(node); i++) {
                singleAddChild(new Node(tree, tree.getChild(node, i)));
            }
            // System.err.println("Built initial tree");
            // System.exit(-1);
        }

        // public Node(Node node, int partition)
        // {
        // leftParent = rightParent = null;
        // leftChild = rightChild = null;
        // leftPartition = node.leftPartition;
        // rightPartition = node.rightPartition;
        // //leftPartition = rightPartition = null;
        // heightParameter = node.heightParameter;
        // number = node.number;
        // taxon = node.taxon;
        // bifurcation = node.bifurcation;
        // // System.err.println("Examinging "+number);
        // Node left;
        // if( node.leftChild != null ) {
        // Node left = node.getChild(0,partition);
        // }
        // Node right;
        // if( node.rightChild != null ) {
        // Node right = node.getChild(1,partition);
        // }
        //
        // if( left != null ) {
        // // System.err.println("Adding "+number+"->"+left.number);
        // singleAddChild(new Node(left,partition));
        // }
        // }
        //
        // if( right != null ) {
        // // System.err.println("Adding "+number+"->"+right.number);
        // singleAddChild(new Node(right,partition));
        // }
        // }
        // }

        public Node(Node inode, int partition) { //, ArrayList<Node> nodes) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            Node node = inode;
            while (node.isBifurcationDoublyLinked()) {
                node = node.leftChild.leftChild;
                // System.err.println("Does this do anything?");
            }
            heightParameter = node.heightParameter;
            number = node.number;
            taxon = node.taxon;
            bifurcation = true;

            mirrorNode = node;

            if (node.isExternal())
//                nodes.add(this);
                return;
            else {
                Node left, right;
                left = node.getChild(0, partition);
                right = node.getChild(1, partition);
                if (left != null || right != null) {
                    if (left != null) {
                        singleAddChild(new Node(left, partition));
                    }
                    if (right != null) {
                        singleAddChild(new Node(right, partition));
                    }
//                  nodes.add(this);
                }
            }
        }

        public void setPartitionRecursively(int partition) {
            boolean onLeft = MathUtils.nextBoolean();
            if (leftChild != null) {
                // if( leftPartition != null )
                // leftPartition.set(partition);
                if (partitioning != null && onLeft)
                    partitioning.setParameterValue(partition, 1);
                leftChild.setPartitionRecursively(partition);
            }
            if (rightChild != null) {
                // if( leftPartition != null )
                // rightPartition.set(partition);
                if (partitioning != null && !onLeft)
                    partitioning.setParameterValue(partition, 1);
                rightChild.setPartitionRecursively(partition);
            }

        }

        public void stripOutDeadEnds() {
            if (leftChild != null)
                leftChild.stripOutDeadEnds();
            if (rightChild != null && rightChild != leftChild)
                rightChild.stripOutDeadEnds();
            if (taxon == null && leftChild == null && rightChild == null)
                leftParent.doubleRemoveChild(this);
        }

        public Node stripOutSingleChildNodes(Node cRoot) {
            // Node rtn = cRoot;
            int childCount = getChildCount();
            if (childCount == 0) {
                return cRoot;
            }
            if (childCount == 2) {
                if (hasEqualChildren()) {
                    return leftChild.stripOutSingleChildNodes(leftChild);
                }
                leftChild.stripOutSingleChildNodes(cRoot);
                rightChild.stripOutSingleChildNodes(cRoot);
                return cRoot;
            }
            if (isRoot()) {
                if (leftChild != null) {
                    leftChild.leftParent = leftChild.rightParent = null;
                    return leftChild.stripOutSingleChildNodes(leftChild);
                } else {
                    rightChild.leftParent = rightChild.rightParent = null;
                    return rightChild.stripOutSingleChildNodes(rightChild);
                }
            }
            // System.err.println("Unlinking "+number);
            Node parent = leftParent;
            Node child = leftChild;
            if (child == null)
                child = rightChild;
            parent.doubleRemoveChild(this);
            doubleRemoveChild(child);
            parent.singleAddChild(child);
            return child.stripOutSingleChildNodes(cRoot);
        }

        public final void setupHeightBounds() {
            heightParameter.addBounds(new NodeHeightBounds(heightParameter));
        }


        public final void createRateParameter(int numberPartitions) {
            if (rateParameter == null) {

                double[] startingRateValues = new double[numberPartitions];

                for (int i = 0; i < startingRateValues.length; i++) {
                    startingRateValues[i] = 1.0;
                }

                rateParameter = new Parameter.Default(startingRateValues);

                if (isRoot()) {
                    rateParameter.setId("root.rate");
                } else if (isExternal()) {
                    rateParameter.setId(getTaxonId(getNumber()) + ".rate");
                } else {
                    rateParameter.setId("node" + getNumber() + ".rate");
                }
                rateParameter.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, 0, startingRateValues.length));
                addVariable(rateParameter);
            }
        }

        public final void createTraitParameter() {
            if (traitParameter == null) {
                traitParameter = new Parameter.Default(1.0);
                if (isRoot()) {
                    traitParameter.setId("root.trait");
                } else if (isExternal()) {
                    traitParameter.setId(getTaxonId(getNumber()) + ".trait");
                } else {
                    traitParameter.setId("node" + getNumber() + ".trait");
                }
                rateParameter.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

                addVariable(traitParameter);
            }
        }

        public final double getHeight() {
            return heightParameter.getParameterValue(0);
        }

        public final double getRate(int partition) {


            return rateParameter.getParameterValue(partition);
        }

        public final double getTrait() {
            return traitParameter.getParameterValue(0);
        }

        public final void setHeight(double height) {
            heightParameter.setParameterValue(0, height);
        }

        public final void setRate(double rate) {
            // System.out.println("Rate set for parameter " +
            // rateParameter.getParameterName());
            rateParameter.setParameterValue(0, rate);
        }

        public final void setTrait(double trait) {
            // System.out.println("Trait set for parameter " +
            // traitParameter.getParameterName());
            traitParameter.setParameterValue(0, trait);
        }

        public int getNumber() {
            return number;
        }

        public void setNumber(int n) {
            number = n;
        }

        /**
         * Returns the number of children this node has.
         */
        public final int getChildCount() {
            int n = 0;
            if (leftChild != null)
                n++;
            if (rightChild != null && bifurcation)
                n++;
            return n;
        }

        public final int getParentCount() {
            int n = 0;
            if (leftParent != null)
                n++;
            if (rightParent != null)
                n++;
            return n;
        }

        public Node getChild(int n) {
            if (n == 0)
                return leftChild;
            if (n == 1)
                return rightChild;
            throw new IllegalArgumentException(
                    "ARGModel.Nodes can only have up to 2 children");
        }


        private boolean isBifurcationDoublyLinked() {
            return bifurcation && (leftChild == rightChild)
                    && leftChild != null;
        }

        private boolean recombinantIsLinked(Node parent, int partition) {
            boolean left = leftParent == parent;
            boolean right = rightParent == parent;
            final double partitionSide = partitioning
                    .getParameterValue(partition);
            if (left && partitionSide == 0)
                return true;
            if (right && partitionSide == 1)
                return true;
            return false;

        }

        public int getDescendentTipCount() {
            if (isExternal())
                return 1;
            return leftChild.getDescendentTipCount()
                    + rightChild.getDescendentTipCount();
        }

        public boolean checkForNullRights() {
            if (isExternal())
                return false;
            if (rightChild == null)
                return true;
            else
                return rightChild.checkForNullRights()
                        || leftChild.checkForNullRights();
        }

        private Node findNextTreeNode(Node parent, int partition) { // searches
            // down the
            // ARG for
            // the next
            // bifurcation
            // or tip
            if (isExternal())
                return this;
            Node next = this;
            while (next.isReassortment()) {
                if (recombinantIsLinked(parent, partition))
                    return leftChild.findNextTreeNode(parent, partition);
                else
                    return null;
            }
            // if( leftChild.findNextTreeMode())
            next = leftChild.findNextTreeNode(parent, partition);
            if (next == null)
                next = rightChild.findNextTreeNode(parent, partition);
            if (next == null) // TODO Error check for removal
                throw new IllegalArgumentException("Can't find next tree node.");
            return next;
        }

        private Node getChild(int n, int partition) { // Assuming an acyclic
            // bifurcating tree
            Node child = null;
            if (n == 0) // Handle left side
                child = leftChild;
            if (n == 1) // Handle right side
                child = rightChild;
            if (child.isExternal())
                return child;
            if (child.isBifurcationDoublyLinked()) {
                // System.err.println("Passing double from
                // "+number+"->"+child.leftChild.number);
                return child.leftChild.getChild(0, partition);
            }
            if (child.isReassortment()) {
                if (child.recombinantIsLinked(this, partition))
                    return child.getChild(0, partition);
                else
                    return null;
            }
            if (child.leftChild.isReassortment()
                    && child.rightChild.isReassortment()) {
                if (child.leftChild.recombinantIsLinked(child, partition))
                    return child;
                if (child.rightChild.recombinantIsLinked(child, partition))
                    return child;
                return null;
            }
            return child;
        }

        public Node getParent(int n) {
            if (n == 0)
                return leftParent;
            if (n == 1)
                return rightParent;
            throw new IllegalArgumentException(
                    "ARGModel.Nodes can only have 2 parents");
        }

        public boolean hasChild(Node node) {
            return (leftChild == node || rightChild == node);
        }

        /**
         * add new child node
         *
         * @param node new child node
         */
        public void singleAddChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException(
                        "ARGModel.Nodes can only have 2 children");
            }
            // if( node.leftParent == null )
            if (node.leftParent == null)
                node.leftParent = this;
            if (node.rightParent == null)
                node.rightParent = this;
        }

        public void singleAddChildWithOneParent(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException("ARGModel.Node " + number
                        + " can only have 2 children");
            }
            // if( node.leftParent == null )
            if (node.leftParent == null) {
                node.leftParent = this;
            } else if (node.rightParent == null) {
                node.rightParent = this;
            } else {
                throw new IllegalArgumentException(
                        "ARGModel.Nodes can only have 2 parents");
            }
        }

        public void doubleAddChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            }
            if (rightChild == null) {
                rightChild = node;
            }// else {
            // throw new IllegalArgumentException("ARGModel.Nodes can only have
            // 2 children");
            // }
            // if( node.leftParent == null )
            if (node.leftParent == null)
                node.leftParent = this;
            if (node.rightParent == null)
                node.rightParent = this;
        }

        public void doubleAddChildWithOneParent(Node node) {
            if (leftChild == null) {
                leftChild = node;
            }
            if (rightChild == null) {
                rightChild = node;
            }// else {
            // throw new IllegalArgumentException("ARGModel.Nodes can only have
            // 2 children");
            // }
            if (node.leftParent == null) {
                node.leftParent = this;
            } else if (node.rightParent == null) {
                node.rightParent = this;
            } else {
                throw new IllegalArgumentException(
                        "ARGModel.Nodes can only have 2 parents");
            }
        }

        public void addChildNoParentConnection(Node node) {
            if (leftChild == null)
                leftChild = node;
            else if (rightChild == null)
                rightChild = node;
            else
                throw new IllegalArgumentException(
                        "Nodes can only have 2 children.");
        }

        // public void addChild()

        // public void addChildRecombinant(Node node, BitSet partition) {

        public void addChildRecombinant(Node node, Parameter partition) {
            // if( leftChild == null && rightChild == null ) {
            // leftChild = rightChild = node;
            // } else
            // if( leftChild == null && rightChild == null ) {
            // System.err.println("yep");
            // //System.exit(-1);
            // leftChild = rightChild = node;
            // }
            if (leftChild == null)
                leftChild = node;
            // leftPartition = partition;
            if (rightChild == null)
                rightChild = node;
            // rightPartition = partition;
            // } else {
            // throw new IllegalArgumentException("Nodes can only have 2
            // children.");
            // }
            // node.parent = null;
            if (node.leftParent == null) {
                node.leftParent = this;
                // node.leftPartition = partition;
                node.partitioning = partition;
            } else if (node.rightParent == null) {
                node.rightParent = this;
                // node.rightPartition = partition;
                node.partitioning = partition;
            } else {
                throw new IllegalArgumentException(
                        "Recombinant nodes can only have 2 parents.");
            }
        }

        /**
         * remove child
         *
         * @param node child to be removed
         */
        public Node doubleRemoveChild(Node node) {
            boolean found = false;
            if (leftChild == node) {
                leftChild = null;
                // leftPartition = null;
                found = true;
            }
            if (rightChild == node) {
                rightChild = null;
                // rightPartition = null;
                found = true;
            }
            if (!found)
                throw new IllegalArgumentException("Unknown child node");
            if (node.leftParent == this)
                node.leftParent = null;
            // node.leftPartition = null;
            if (node.rightParent == this)
                node.rightParent = null;
            // node.rightPartition = null;
            return node;
        }

        public Node singleRemoveChild(Node node) {
            // boolean found = false;
            if (leftChild == node) {
                leftChild = null;
                // leftPartition = null;
                // found = true;
            } else if (rightChild == node) {
                rightChild = null;
                // rightPartition = null;
                // found = true;
            }
            // if( !found )
            else
                throw new IllegalArgumentException("Unknown child node");
            if (node.bifurcation) {
                node.leftParent = node.rightParent = null;
                return null;
            }
            if (node.leftParent == this)
                node.leftParent = null;
                // node.leftPartition = null;
            else if (node.rightParent == this)
                node.rightParent = null;
            // node.rightPartition = null;
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
                throw new IllegalArgumentException(
                        "TreeModel.Nodes can only have 2 children");
            }
            if (node.leftParent == this)
                node.leftParent = null;
            if (node.rightParent == this)
                node.rightParent = null;
            return node;
        }

        public boolean hasChildren() {
            return (leftChild != null || rightChild != null);
        }

        public boolean isExternal() {
            return !hasChildren();
        }

        public boolean isRoot() {
            return (leftParent == null && rightParent == null);
        }

        public boolean hasEqualChildren() {
            return (leftChild == rightChild);
        }

        public boolean isBifurcation() {
            return bifurcation;
        }

        // public boolean isReassortment() { return hasChildren() && (leftChild
        // == rightChild); }
        public boolean isReassortment() {
            return !bifurcation;
        }

        private String toExtendedNewick() {

            if (isExternal())
                return taxon.getId();

            if (isBifurcation()) {
                String left = leftChild.toExtendedNewick();
                String right = rightChild.toExtendedNewick();
                if (left.compareTo(right) < 0)
                    return "(" + left + "," + right + ")";
                else
                    return "(" + right + "," + left + ")";
            }
            // must be a reassortment node
            return "<" + leftChild.toExtendedNewick() + ">";
        }

        public String toString() {
            if (taxon == null) {
                return "" + number;
            }
            return "" + number + " (" + taxon.getId() + ")";
        }
    }

    /**
     * This class provides bounds for parameters that represent a node height in
     * this tree model.
     */
    private class NodeHeightBounds implements Bounds<Double> {

        public NodeHeightBounds(Parameter parameter) {
            nodeHeightParameter = parameter;
        }

        public Double getUpperLimit(int i) {
            // I think only upper bounds are of concern with linked subtrees
            // because everything below has only one parameter
            // TODO -- check this!
            Node node = getNodeOfParameter(nodeHeightParameter);
            // Returns the first node in nodes[] with this height parameter

            if (node.isRoot()) {
                return Double.POSITIVE_INFINITY;
                // return 10.0;
            } else {
                if (node.leftParent == null) {
                    System.err.println("leftParent of " + node.number
                            + " is null");
                }
                if (node.rightParent == null) {
                    System.err.println("rightParent of " + node.number
                            + " is null");
                }
                return Math.min(node.leftParent.getHeight(), node.rightParent
                        .getHeight());
            }
        }

        public Double getLowerLimit(int i) {

            Node node = getNodeOfParameter(nodeHeightParameter);
            // System.err.println("Is node recombinant?
            // "+node.isReassortment());
            if (node.isExternal()) {
                return 0.0;
            } else {
                if (node.leftChild == null)
                    System.err.println("Node " + node.number
                            + " has null leftChild");
                if (node.rightChild == null)
                    System.err.println("Node " + node.number
                            + " has null rightChild");
                // System.err.println(node.number+" "+(node.leftChild==null)+"
                // "+(node.rightChild==null));
                return Math.max(node.leftChild.getHeight(), node.rightChild
                        .getHeight());
            }
        }

        public int getBoundsDimension() {
            return 1;
        }

        private Parameter nodeHeightParameter = null;
    }

    public double getNodeRate(NodeRef node) {
        if (true)
            throw new RuntimeException("This should not be called");

        return 0;
    }

    ///////////////////////////////////////////////////////////////////////
    //PARSER
    ///////////////////////////////////////////////////////////////////////

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_MODEL;
        }

        public String[] getParserNames() {
            return new String[]{
                    getParserName(), "argModel"
            };
        }

        /**
         * @return a tree object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Tree tree = (Tree) xo.getChild(Tree.class);
            ARGModel treeModel = new ARGModel(tree);

            Logger.getLogger("dr.evomodel").info("Creating the tree model, '" + xo.getId() + "'");

            if (xo.hasAttribute(PARTITION_TYPE)) {
                treeModel.partitionType = xo.getStringAttribute(PARTITION_TYPE);

                if (!treeModel.partitionType.equals(REASSORTMENT_PARTITION) &&
                        !treeModel.partitionType.equals(RECOMBINATION_PARTITION)) {
                    throw new XMLParseException("Must use either correct partition type");
                }
            }

            int numberPartitions = 1;

            if (xo.hasAttribute(NUM_PARTITIONS)) {
                numberPartitions = xo.getIntegerAttribute(NUM_PARTITIONS);
            }


            Logger.getLogger("dr.evomodel").info(
                    xo.getId() + " has partition type: " + treeModel.partitionType);

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof XMLObject) {

                    XMLObject cxo = (XMLObject) xo.getChild(i);

                    if (cxo.getName().equals(ROOT_HEIGHT)) {

                        ParameterParser.replaceParameter(cxo, treeModel
                                .getRootHeightParameter());

                    } else if (cxo.getName().equals(LEAF_HEIGHT)) {

                        String taxonName;
                        if (cxo.hasAttribute(TAXON)) {
                            taxonName = cxo.getStringAttribute(TAXON);
                        } else {
                            throw new XMLParseException(
                                    "taxa element missing from leafHeight element in treeModel element");
                        }

                        int index = treeModel.getTaxonIndex(taxonName);
                        if (index == -1) {
                            throw new XMLParseException(
                                    "taxon "
                                            + taxonName
                                            + " not found for leafHeight element in treeModel element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);
                        ParameterParser.replaceParameter(cxo, treeModel
                                .getLeafHeightParameter(node));

                    } else if (cxo.getName().equals(NODE_HEIGHTS)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo
                                    .getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException(
                                    "one or more of root, internal or leaf nodes must be " +
                                            "selected for the nodeHeights element");
                        }

                        ParameterParser.replaceParameter(cxo, treeModel
                                .createNodeHeightsParameter(rootNode,
                                internalNodes, leafNodes));

                    } else if (cxo.getName().equals(NODE_RATES)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo
                                    .getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        // if (rootNode) {
                        // throw new XMLParseException("root node does not have
                        // a rate parameter");
                        // }

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException(
                                    "one or more of root, internal or leaf nodes must be selected for the nodeRates element");
                        }

                        ParameterParser.replaceParameter(cxo, treeModel
                                .createNodeRatesParameter(rootNode,
                                internalNodes, leafNodes, numberPartitions));

                    } else if (cxo.getName().equals(NODE_TRAITS)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo
                                    .getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException(
                                    "one or more of root, internal or leaf nodes must be selected for the nodeTraits element");
                        }

                        ParameterParser.replaceParameter(cxo, treeModel
                                .createNodeTraitsParameter(rootNode,
                                internalNodes, leafNodes));

                    } else {
                        throw new XMLParseException("illegal child element in "
                                + getParserName() + ": " + cxo.getName());
                    }

                } else if (xo.getChild(i) instanceof Tree) {
                    // do nothing - already handled
                } else {
                    throw new XMLParseException("illegal child element in  "
                            + getParserName() + ": " + xo.getChildName(i) + " "
                            + xo.getChild(i));
                }
            }


            treeModel.setupHeightBounds();

            Logger.getLogger("dr.evomodel").info(
                    "  initial tree topology = "
                            + TreeUtils.uniqueNewick(treeModel, treeModel
                            .getRoot()));
            return treeModel;
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents a model of the tree. The tree model includes and attributes of the nodes "
                    + "including the age (or <i>height</i>) and the rate of evolution at each node in the tree.";
        }

        public String getExample() {
            return "<!-- the tree model as special sockets for attaching parameters to various aspects of the tree     -->\n"
                    + "<!-- The treeModel below shows the standard setup with a parameter associated with the root height -->\n"
                    + "<!-- a parameter associated with the internal node heights (minus the root height) and             -->\n"
                    + "<!-- a parameter associates with all the internal node heights                                     -->\n"
                    + "<!-- Notice that these parameters are overlapping                                                  -->\n"
                    + "<!-- The parameters are subsequently used in operators to propose changes to the tree node heights -->\n"
                    + "<treeModel id=\"treeModel1\">\n"
                    + "	<tree idref=\"startingTree\"/>\n"
                    + "	<rootHeight>\n"
                    + "		<parameter id=\"treeModel1.rootHeight\"/>\n"
                    + "	</rootHeight>\n"
                    + "	<nodeHeights internalNodes=\"true\" rootNode=\"false\">\n"
                    + "		<parameter id=\"treeModel1.internalNodeHeights\"/>\n"
                    + "	</nodeHeights>\n"
                    + "	<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n"
                    + "		<parameter id=\"treeModel1.allInternalNodeHeights\"/>\n"
                    + "	</nodeHeights>\n" + "</treeModel>";

        }

        public Class getReturnType() {
            return ARGModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final String[] partitionFormats = {REASSORTMENT_PARTITION, RECOMBINATION_PARTITION};

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(PARTITION_TYPE, "Describes the partition structure of the model",
                        partitionFormats, true),
                new ElementRule(Tree.class),
                new ElementRule(
                        ROOT_HEIGHT,
                        Parameter.class,
                        "A parameter definition with id only (cannot be a reference!)",
                        false),
                new ElementRule(
                        NODE_HEIGHTS,
                        new XMLSyntaxRule[]{
                                AttributeRule
                                        .newBooleanRule(ROOT_NODE, true,
                                        "If true the root height is included in the parameter"),
                                AttributeRule
                                        .newBooleanRule(
                                        INTERNAL_NODES,
                                        true,
                                        "If true the internal node heights (minus the root) are included in the parameter"),
                                new ElementRule(Parameter.class,
                                        "A parameter definition with id only (cannot be a reference!)")},
                        1, Integer.MAX_VALUE)};

        public Parameter oldGetParameter(XMLObject xo) throws XMLParseException {
            //		public Parameter getParameter(XMLObject xo) throws XMLParseException {


            int paramCount = 0;
            Parameter param = null;
            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Parameter) {
                    param = (Parameter) xo.getChild(i);
                    paramCount += 1;
                }
            }

            if (paramCount == 0) {
                throw new XMLParseException(
                        "no parameter element in treeModel " + xo.getName()
                                + " element");
            } else if (paramCount > 1) {
                throw new XMLParseException(
                        "More than one parameter element in treeModel "
                                + xo.getName() + " element");
            }

            return param;
        }

        // todo check to make sure that Andrew's static routine works with this old code
        public void oldReplaceParameter(XMLObject xo, Parameter newParam)
//		public void replaceParameter(XMLObject xo, Parameter newParam)

                throws XMLParseException {

            for (int i = 0; i < xo.getChildCount(); i++) {

                if (xo.getChild(i) instanceof Parameter) {

                    XMLObject rxo = null;
                    Object obj = xo.getRawChild(i);

                    if (obj instanceof Reference) {
                        rxo = ((Reference) obj).getReferenceObject();
                    } else if (obj instanceof XMLObject) {
                        rxo = (XMLObject) obj;
                    } else {
                        throw new XMLParseException(
                                "object reference not available");
                    }

                    if (rxo.getChildCount() > 0) {
                        throw new XMLParseException(
                                "No child elements allowed in parameter element.");
                    }

                    if (rxo.hasAttribute(XMLParser.IDREF)) {
                        throw new XMLParseException("References to "
                                + xo.getName()
                                + " parameters are not allowed in treeModel.");
                    }

                    if (rxo.hasAttribute(ParameterParser.VALUE)) {
                        throw new XMLParseException("Parameters in "
                                + xo.getName()
                                + " have values set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.UPPER)) {
                        throw new XMLParseException("Parameters in "
                                + xo.getName()
                                + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.LOWER)) {
                        throw new XMLParseException("Parameters in "
                                + xo.getName()
                                + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(XMLParser.ID)) {

                        newParam.setId(rxo.getStringAttribute(XMLParser.ID));
                    }

                    rxo.setNativeObject(newParam);

                    return;
                }
            }
        }
    };


}
