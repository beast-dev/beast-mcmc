package dr.evomodel.bigFastTree;

/*
 * SubtreeLeapOperator.java
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


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.operators.AbstractAdaptableTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.bigFastTree.CladeAwareSubtreeLeapOperatorParser;
import dr.inference.distribution.CauchyDistribution;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.TruncatedNormalDistribution;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements a clade aware Subtree Leap move.
 *
 * This move picks a node at random (except for the root) and then moves the parent to any location
 * that is a certain patristic distance from its starting point without breaking a defined clade structure .
 *
 * It is always possible for the node in the root clade to move up (potentially becoming the root) but the destination can't
 * be younger than the original node. All possible destinations are collected and then picked amongst
 * uniformly.
 *
 * @author Andrew Rambaut
 * @author Luiz Max Carvalho
 * @author Mathieu Fourment
 * @author JT McCrone
 * @version $Id$
 */
public class CladeAwareSubtreeLeap extends AbstractAdaptableTreeOperator {


    public enum DistanceKernelType {
        NORMAL("normal") {
            @Override
            double getDelta(double size, double bound) {
                Distribution D;
                if(bound == Double.POSITIVE_INFINITY){
                    D = new NormalDistribution(0,size);
                }else{
                    D = new TruncatedNormalDistribution(0, size, 0, bound);
                }
                double u = MathUtils.nextDouble();
                return  Math.abs(D.quantile(u));
            }

            @Override
            double getPdf(double size, double bound, double delta) {
                Distribution D;
                if(bound == Double.POSITIVE_INFINITY){
                      D = new NormalDistribution(0,size);
                }else{
                     D = new TruncatedNormalDistribution(0, size, 0, bound);
                }
                return D.pdf(delta);
            }
        },
        CAUCHY("Cauchy") {
            @Override
            double getDelta(double size,double bound) {
                Distribution distK = new CauchyDistribution(0, size);
                double u = MathUtils.nextDouble();
                return Math.abs(distK.quantile(u));
            }
            @Override
            double getPdf(double size, double bound, double delta) {
                throw new RuntimeException("Cauchy not implemented for CladeAware SubtreeLeap");
//                return 0;
            }
        };

        DistanceKernelType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;

        abstract double getDelta(double size,double bound);

        abstract double getPdf(double size, double bound, double delta);
    }

    private double size;

    private final TreeModel tree;
    private final DistanceKernelType distanceKernel;
    private final boolean slideOnly;
    private final CladeNodeModel cladeModel;
    private final boolean DEBUG=false;

    /**
     * Constructor
     * @param cladeModel   the cladeModel
     * @param weight the weight
     */
    public CladeAwareSubtreeLeap(CladeNodeModel cladeModel, double weight)  {
        this(cladeModel, weight, 1.0, DistanceKernelType.NORMAL,AdaptationMode.DEFAULT,0.225, false);
    }
    /**
     * Constructor
     * @param weight the weight
     */
    public CladeAwareSubtreeLeap(CladeNodeModel cladeModel, double weight,boolean slideOnly) {
        this(cladeModel, weight, 1.0, DistanceKernelType.NORMAL, AdaptationMode.DEFAULT,0.225,slideOnly);
    }
    /**
     * Constructor
     *
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param distanceKernel the distribution from which to draw the patristic distance
     * @param slideOnly if true, only slide up and down the tree, never across (mimics SubtreeSlide)
     */
    public CladeAwareSubtreeLeap(CladeNodeModel cladeModel, double weight, double size, DistanceKernelType distanceKernel, AdaptationMode mode, double targetAcceptance, boolean slideOnly) {
        super(mode, targetAcceptance);

        this.cladeModel = cladeModel;
        this.tree = cladeModel.getTreeModel();
        setWeight(weight);
        this.size = size;
        this.distanceKernel = distanceKernel;
        this.slideOnly = slideOnly;

    }


    public double getMaxDistance(NodeRef node, NodeRef parent, NodeRef sibling ) {
        CladeRef clade = cladeModel.getClade(parent);


        double height = tree.getNodeHeight(node);

        double distanceBelowParent = searchDown(parent, sibling, clade, height);

        NodeRef node1 = parent;

        // walk up to clade root
        double distanceAboveParent=0;
        double distanceUp = 0;
        double currentMax = distanceBelowParent;

        boolean done = false;
        while (!done) {
            NodeRef parent1 = tree.getParent(node1);
            if (parent1 != null && cladeModel.getClade(parent1) == clade) {
                 distanceUp = tree.getNodeHeight(parent1)-tree.getNodeHeight(parent);
                 currentMax = Math.max(distanceUp, currentMax);
                    if (!slideOnly) { // if we are not just sliding up or down...
                        double distanceDown;
                        // We haven't reached the height above the original height so go down
                        // the sibling subtree to look for other possible destinations
                        NodeRef sibling1 = getOtherChild(tree, parent1, node1);

                        if (tree.getNodeHeight(sibling1) > height) {
                            distanceDown = searchDown(parent1, sibling1, clade, height);
                        }else{
                            distanceDown = tree.getNodeHeight(parent1) - height;
                        }

                        double subtreeDistance = distanceUp+distanceDown;
                        currentMax =  Math.max(subtreeDistance, currentMax);
                    }

                node1 = parent1;
            } else {
                // node1 is the root - add it as a destination and stop loop
                if (parent1 == null) {
                    throw new RuntimeException("Should never hit the root in this traversal");
//                    done = true;
                } else if (cladeModel.getClade(parent1) != clade) {
                    assert cladeModel.getClade(parent1)== cladeModel.getParent(clade);
                    distanceUp = tree.getNodeHeight(parent1) - tree.getNodeHeight(parent);
                    currentMax = Math.max(distanceUp, currentMax);
                    done = true;
                }
            }
        }
        return currentMax;
    }

    public double searchDown(NodeRef parent, NodeRef node,CladeRef clade, double height) {
        if(tree.getNodeHeight(parent)<height || cladeModel.getClade(parent)!=clade){
            return 0;
        }
        if(tree.getNodeHeight(node)<height){
            return tree.getNodeHeight(parent)-height;
        }
        if (cladeModel.getClade(node) != clade) {
            return tree.getBranchLength(node);
        }
        double distance = tree.getBranchLength(node);
        double[] childDistances = {0,0};
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            childDistances[i] = searchDown(node, child,clade, height);
        }
        return distance + Math.max(childDistances[0],childDistances[1]);
    }



    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq;

        final NodeRef root = tree.getRoot();

        NodeRef node;

            // Pick a node (but not the root)
            do {
                // choose a random node avoiding root
                node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            } while (node == root);

        // get its parent - this is the node we will prune/graft
        final NodeRef parent = tree.getParent(node);

        final CladeRef clade = cladeModel.getClade(parent);

        // get the node's sibling
        final NodeRef sibling = getOtherChild(tree, parent, node);



        double bound =  cladeModel.getParent(clade) == null ? Double.POSITIVE_INFINITY : getMaxDistance(node, parent,sibling);

        double delta = distanceKernel.getDelta(size,bound);

//        while (delta > bound) {
//                delta = distanceKernel.getDelta(size);
//        }


        final Map<NodeRef, Double> destinations = getDestinations(node, parent, sibling, delta, slideOnly);
        final List<NodeRef> destinationNodes = new ArrayList<NodeRef>(destinations.keySet());

        NodeRef[] treenodes = tree.getNodes();
        if(DEBUG) {
            for (NodeRef dn : destinationNodes) {
                boolean found = false;
                for (NodeRef tn : treenodes) {
                    if (tn == dn) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    System.out.println("¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡   destination not in tree  !!!!!!!!!!!!!!!!!!!!!!!!!");
                    System.out.println(dn);
                }
            }
        }
        // pick uniformly from this list
        int r = MathUtils.nextInt(destinations.size());

        double forwardProbability = distanceKernel.getPdf(size,bound,delta) / (destinations.size());

        final NodeRef j = destinationNodes.get(r);
        final double newHeight = destinations.get(j);
        //Make the edits
        doOperation(node,parent, sibling, j, clade, newHeight);

        NodeRef newSibling = getOtherChild(tree, parent, node);

        double reverseBound = cladeModel.getParent(clade) == null ? Double.POSITIVE_INFINITY : getMaxDistance(node, parent,newSibling);

        if(reverseBound<delta){
            throw new RuntimeException("reverse Jump not possible!");
        }

        final Map<NodeRef, Double> reverseDestinations = getDestinations(node, parent,newSibling , delta, slideOnly);
        double reverseProbability = distanceKernel.getPdf(size,reverseBound,delta) / (reverseDestinations.size());

        // hastings ratio = reverse Prob / forward Prob
        logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        return logq;
    }

    /**
     * Do a subtree leap move with specified node, j and height
     *
     * for testing purposes
     *
     */
    public void doOperation(NodeRef node , NodeRef parent, NodeRef sibling, NodeRef j, CladeRef clade,double newHeight) {

        // and its grand parent
        final NodeRef grandParent = tree.getParent(parent);

        final NodeRef jParent = tree.getParent(j);

        if (jParent != null && newHeight > tree.getNodeHeight(jParent)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(j)) {
            throw new IllegalArgumentException("height error");
        }

        if (DEBUG) {
            System.out.println("#############--------- Operation -----------############\"");
            System.out.println("Operating on: " + node);
            System.out.println("Parent: "+ parent);
            System.out.println("Grandparent " + grandParent);
            System.out.println("Sibling " + sibling);
            System.out.println("Destination " + j);
            System.out.println("Destination parent " + jParent);

            System.out.println("New Height " + newHeight);
            System.out.println(tree.toString());

            int[] postOrder = new int[tree.getNodeCount()];
            TreeUtils.postOrderTraversalList(tree,postOrder);
            System.out.print("The postOrder nodes: [ ");
            for (int k : postOrder) {
                System.out.print(" " + k + " ");
            }
            System.out.println("]");

            System.out.print("The tip id:number nodes: [ ");
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                System.out.print(" " + tree.getNodeTaxon(tree.getExternalNode(i)).getId() +":"+ tree.getExternalNode(i).getNumber() +" ");
            }
            System.out.println("]");



        }


        tree.beginTreeEdit();

        // j is the node below where we are going

        if (j == parent || jParent == parent) {
            // the subtree is not actually moving but the height will change
        } else {
            if (grandParent == null) {
                // if the parent of the original node is the root then the sibling becomes
                // the root.
                tree.removeChild(parent, sibling);
                tree.setRoot(sibling);
                assert cladeModel.getParent(clade)==null;
                cladeModel.setRootNode(clade,sibling);

            } else {
                // remove the parent of node by connecting its sibling to its grandparent.
                tree.removeChild(parent, sibling);
                tree.removeChild(grandParent, parent);
                tree.addChild(grandParent, sibling);

            }

            if (jParent == null) {
                // adding the node to the root of the tree
                tree.addChild(parent, j);
                tree.setRoot(parent);
                cladeModel.setRootNode(clade,parent);
                assert cladeModel.getParent(clade)==null;

                if(DEBUG){
                    System.out.println("##############-------supper seed root!-------------###########3");
                    System.out.println("Destination node " + j +" added as child to node:  " + parent);
                    System.out.println("new destination parent = " + tree.getParent(j));

                    System.out.print("parent's children : [ ");
                    for (int i = 0; i < tree.getChildCount(parent); i++) {
                        System.out.print(" " + tree.getChild(parent, i) + " ");
                    }
                    System.out.println("]");

                }
            } else {
                // remove destination edge j from its parent
                tree.removeChild(jParent, j);

                // add destination edge to the parent of node
                tree.addChild(parent, j);

                // and add the parent of i as a child of the former parent of j
                tree.addChild(jParent, parent);

                if (parent == cladeModel.getRootNode(clade)) {
                    cladeModel.setRootNode(clade,sibling);
                }
                // should never be called if the if above is called
                if(cladeModel.getClade(jParent)!=clade){
                    assert cladeModel.getClade(jParent) == cladeModel.getParent(clade);
                    assert cladeModel.getRootNode(clade) == j;
                    cladeModel.setRootNode(clade,parent);

                }
            }

            if (DEBUG) {
                System.out.println("#############--------- after move before edit done -----------############");
                System.out.println("parent : " + parent);
                System.out.print("parent's children : [ ");
                for (int i = 0; i < tree.getChildCount(parent); i++) {
                    System.out.print(" " + tree.getChild(parent, i) + " ");
                }
                System.out.println("]");

                System.out.println(tree.getChild(parent,1)==j);
                System.out.println(tree.getParent(j)==parent);
                System.out.println(j==tree.getNode(j.getNumber()));

                System.out.println("root: "+ tree.getRoot());

            }

        }


        tree.endTreeEdit();
        assert tree==cladeModel.getTreeModel();

        if (DEBUG) {
            System.out.println("#############--------- edit done -----------############");

            System.out.println("parent : " + parent);
            System.out.println("root: "+ tree.getRoot());

            for (int i = 0; i < tree.getChildCount(parent); i++) {
                System.out.print(" " + tree.getChild(parent, i) + " ");
            }
            System.out.println("]");

            System.out.println(tree.getChild(parent,1)==j);
            System.out.println(tree.getParent(j)==parent);

            System.out.println(j==tree.getNode(j.getNumber()));






            for (NodeRef n :tree.getNodes()  ) {
                if (tree.getParent(n) == null) {
                    if (tree.getRoot() !=  n) {
                        System.out.println("Oh no! lost node: " + n.getNumber());
                    }
                }
            }
        }

        tree.setNodeHeight(parent, newHeight);

        if (DEBUG) {
            System.out.println("#############--------- height set -----------############");

            System.out.println(tree.toString());

            int[] postOrder = new int[tree.getNodeCount()];
            TreeUtils.postOrderTraversalList(tree,postOrder);
            System.out.print("The postOrder nodes: [ ");
            for (int k : postOrder) {
                System.out.print(" " + k + " ");
            }
            System.out.println("]");
        }



        if (tree.getParent(parent) != null && newHeight > tree.getNodeHeight(tree.getParent(parent))) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(node)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(getOtherChild(tree, parent, node))) {
            throw new IllegalArgumentException("height error");
        }

    }

    public Map<NodeRef, Double> getDestinations(NodeRef node, NodeRef parent, NodeRef sibling, double delta, boolean slideOnly) {
        final CladeRef clade = cladeModel.getClade(parent);

        final Map<NodeRef, Double> destinations = new LinkedHashMap<NodeRef, Double>();
        // get the parent's height
        final double height = tree.getNodeHeight(parent);

        final double heightBelow = height - delta;

        if (heightBelow > tree.getNodeHeight(node)) {
            // the destination height below the parent is compatible with the node
            // see if there are any destinations on the sibling's branch
            final List<NodeRef> edges = new ArrayList<NodeRef>();

            getIntersectingEdges(sibling, heightBelow, edges, clade);

            // add the intersecting edges and the height
            for (NodeRef n : edges) {
                destinations.put(n, heightBelow);
            }
        }

        final double heightAbove = height + delta;

        NodeRef node1 = parent;

        // walk up to root
        boolean done = false;
        while (!done) {
            NodeRef parent1 = tree.getParent(node1);

            if (parent1 != null && cladeModel.getClade(parent1)==clade) {
                final double height1 = tree.getNodeHeight(parent1);
                if (height1 < heightAbove) {
                    if (!slideOnly) { // if we are not just sliding up or down...

                        // We haven't reached the height above the original height so go down
                        // the sibling subtree to look for other possible destinations
                        NodeRef sibling1 = getOtherChild(tree, parent1, node1);

                        double heightBelow1 = height1 - (heightAbove - height1);

                        if (heightBelow1 > tree.getNodeHeight(node)) {

                            final List<NodeRef> edges = new ArrayList<NodeRef>();

                            getIntersectingEdges(sibling1, heightBelow1, edges, clade);

                            // add the intersecting edges and the height
                            for (NodeRef n : edges) {
                                destinations.put(n, heightBelow1);
                            }
                        }
                    }
                } else {
                    // add the current node as a destination
                    destinations.put(node1, heightAbove);
                    done = true;
                }

                node1 = parent1;
            } else {
                // node1 is the root - add it as a destination and stop loop
                if(parent1==null){
                    destinations.put(node1, heightAbove);
                    done = true;
                } else if (tree.getNodeHeight(parent1)>heightAbove) {
                    destinations.put(node1, heightAbove);
                    done = true;
                }else{
                    //TODO could get here if parent1 one is outside clade. look for spot there if so;
                    if(destinations.size()==0){
                        throw new RuntimeException("Something went wrong no destinations!");
                    }
                    done = true;
                }
            }
        }

        return destinations;
    }


    private int getIntersectingEdges( NodeRef node, double height, List<NodeRef> edges,CladeRef clade) {

        final NodeRef parent = tree.getParent(node);

        if (tree.getNodeHeight(parent) < height || cladeModel.getClade(parent)!=clade) return 0;

        if (tree.getNodeHeight(node) < height) {
            edges.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
                count += getIntersectingEdges(tree.getChild(node, i), height, edges, clade);
        }
        return count;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        setSize(Math.exp(value));
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(getSize());
    }

    @Override
    public double getRawParameter() {
        return getSize();
    }

    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
            return CladeAwareSubtreeLeapOperatorParser.CLADE_AWARE_SUBTREE_LEAP + "(" + tree.getId() + ")";

    }
}
