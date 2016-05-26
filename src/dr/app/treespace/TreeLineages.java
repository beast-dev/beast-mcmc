/*
 * TreeLineages.java
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

package dr.app.treespace;

import dr.xml.Reference;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;

import java.awt.geom.Rectangle2D;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeLineages {

    private final static String[] LOCATIONS = {
            "Africa",
            "USA",
            "Taiwan",
            "China",
            "Russia",
            "Oceania",
            "Asia",
            "Japan",
            "Mexico",
            "South America",
            "Canada",
            "Europe",
            "Southeast Asia",
            "South Korea"
    };

    private final static String[] AIR_COMMUNITIES = {
            "AC1",
            "AC2",
            "AC3",
            "AC4",
            "AC5",
            "AC6",
            "AC7",
            "AC8",
            "AC9",
            "AC10",
            "AC11",
            "AC12",
            "AC13",
            "AC14"
    };

    Map<String, Integer> locationMap = new HashMap<String, Integer>();
    private static final double THRESHOLD = 0.01;

    public TreeLineages() {
        for (int i = 0; i < AIR_COMMUNITIES.length; i++) {
            locationMap.put(AIR_COMMUNITIES[i], i);
        }
    }

    public List<Lineage> getRootLineages() {
        return rootLineages;
    }

    public Rectangle2D getTreeBounds() {
        return treeBounds;
    }

    public Rectangle2D getAntigenicBounds() {
        return antigenicBounds;
    }

    public void addTree(RootedTree tree) {

        // the offset is the distance to the right most tip - used to align the tips of the tree
        offsetX = 0.0;
        Lineage lineage = new Lineage();
        addNode(tree, tree.getRootNode(), lineage, 0.0);

        lineage.dx = offsetX;

        rootLineages.add(lineage);
    }

    public void setupTrees() {
        for (BitSet key : nodeCounts.keySet()) {
            int count = nodeCounts.get(key);
            if ((((double)count) / rootLineages.size()) > THRESHOLD) {
                majorityNodes.add(key);
            }
        }
        Collections.sort(majorityNodes, new Comparator<BitSet>() {
            public int compare(BitSet bitSet, BitSet bitSet1) {
                return bitSet.cardinality() - bitSet1.cardinality();
            }
        });

        for (Lineage rootLineage : rootLineages) {
        // set the tip count to zero for this traversal
//        currentY = 0.0;
            positionNode(rootLineage);
        }
    }

    public BitSet addNode(RootedTree tree, Node node, Lineage lineage, double cumulativeX) {
        lineage.bits = new BitSet();

        cumulativeX += tree.getLength(node);
        lineage.dx = cumulativeX;

        String location = (String)node.getAttribute("states");
        if (location != null) {
            lineage.state = locationMap.get(location);
        }

        if (!tree.isExternal(node)) {
            List<Node> children = tree.getChildren(node);
            if (children.size() != 2) {
                throw new RuntimeException("Tree is not binary");
            }

            lineage.child1 = new Lineage();
            lineage.child2 = new Lineage();

            lineage.bits.or(addNode(tree, children.get(0), lineage.child1, cumulativeX));
            lineage.bits.or(addNode(tree, children.get(1), lineage.child2, cumulativeX));

            lineage.tipCount = lineage.child1.tipCount + lineage.child2.tipCount;

            if (lineage.child1.tipCount > lineage.child2.tipCount ||
                    (lineage.child1.tipCount == lineage.child2.tipCount &&
                            lineage.child1.tipNumber > lineage.child2.tipNumber)) {
                Lineage tmp = lineage.child1;
                lineage.child1 = lineage.child2;
                lineage.child2 = tmp;
            }

            Integer count = nodeCounts.get(lineage.bits);
            if (count == null) {
                count = 1;
            } else {
                count ++;
            }
            nodeCounts.put(lineage.bits, count);

        } else {
            Integer tipNumber = taxonNumbers.get(tree.getTaxon(node).getName());
            if (tipNumber == null) {
                tipNumber = taxonNumbers.size();
                taxonNumbers.put(tree.getTaxon(node).getName(), tipNumber);
            }

            lineage.tipNumber = tipNumber;
            lineage.tipCount = 1;

            lineage.bits.set(lineage.tipNumber);

            if (cumulativeX > offsetX) {
                offsetX = cumulativeX;
            }
        }

        Object[] ag = (Object[])node.getAttribute("antigenic");
        if (ag != null) {
//            lineage.ax = (Double)ag[0];
//            lineage.ay = (Double)ag[1];
            lineage.ax = lineage.dx;
            lineage.ay = (Double)ag[0];

            antigenicBounds.add(lineage.ax, lineage.ay);
        } else {
            Double ag1 = (Double)node.getAttribute("antigenic1");
            Double ag2 = (Double)node.getAttribute("antigenic2");
            if (ag1 != null) {
                double ad =

                lineage.ay = ag1;
//                lineage.ay = ag2;
                lineage.ax = lineage.dx;

                antigenicBounds.add(lineage.ax, lineage.ay);
            }

        }

        return lineage.bits;
    }

    public double positionNode(Lineage lineage) {
        if (lineage.child1 != null) {
            lineage.dy = positionNode(lineage.child1);
            lineage.dy += positionNode(lineage.child2);

            // the y of this node is the average of the two children
            lineage.dy /= 2.0;

            // now change the children to relative y positions.
//            lineage.child1.dy = lineage.child1.dy - lineage.dy;
//            lineage.child2.dy = lineage.child2.dy - lineage.dy;

            Double location = nodeLocations.get(lineage.bits);
            if (location == null) {
                BitSet bestBits = lineage.bits;

                if (!majorityNodes.contains(bestBits)) {
                    // otherwise find a bigger node...
                    for (BitSet bits : majorityNodes) {
                        BitSet bits1 = (BitSet)bits.clone();
                        bits1.and(lineage.bits);
                        if (bits1.cardinality() == lineage.bits.cardinality()) {
                            bestBits = bits;
                            break;
                        }
                    }

                    location = nodeLocations.get(bestBits);
                    if (location != null) {
                        lineage.dy = location;
                    }
                }

                nodeLocations.put(bestBits, lineage.dy);

            } else {
                lineage.dy = location;
            }

        } else {
            Integer orderedTipNumber = orderedTipNumbers.get(lineage.tipNumber);
            if (orderedTipNumber == null) {
                orderedTipNumber = orderedTipNumbers.size();
                orderedTipNumbers.put(lineage.tipNumber, orderedTipNumber);
            }

            lineage.tipNumber = orderedTipNumber;
//            lineage.dy = LATITUDES[lineage.state];
            lineage.dy = lineage.tipNumber;
//            lineage.dy = currentY;
//            currentY += 1.0;
        }

        treeBounds.add(lineage.dx, lineage.dy);

        return lineage.dy;
    }

    class Lineage {
        double dx = 0;
        double dy = 0;
        double ax = 0;
        double ay = 0;
        Lineage child1 = null;
        Lineage child2 = null;
        int tipNumber = 0;
        int tipCount = 0;
        int state;
        BitSet bits;
    }

    private List<Lineage> rootLineages = new ArrayList<Lineage>();
    private Map<String, Integer> taxonNumbers = new HashMap<String, Integer>();
    private Map<Integer, Integer> orderedTipNumbers = new HashMap<Integer, Integer>();
    private Map<BitSet, Integer> nodeCounts = new HashMap<BitSet, Integer>();
    private Map<BitSet, Double> nodeLocations = new HashMap<BitSet, Double>();
    private List<BitSet> majorityNodes = new ArrayList<BitSet>();

    private double offsetX;

    private Rectangle2D treeBounds = new Rectangle2D.Double();
    private Rectangle2D antigenicBounds = new Rectangle2D.Double();
}
