/*
 * BranchRatePlotter.java
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

package dr.app.tools;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.app.gui.tree.JTreeDisplay;
import dr.app.gui.tree.JTreePanel;
import dr.app.gui.tree.SquareTreePainter;
import dr.stats.DiscreteStatistics;

import javax.swing.*;
import java.awt.*;
import java.awt.print.PrinterJob;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Set;
import java.util.StringTokenizer;

/**
 * Date: Nov 30, 2004
 * Time: 5:29:36 PM
 *
 * @author Alexei Drummond
 *
 * @version $Id: BranchRatePlotter.java,v 1.6 2005/12/08 13:52:46 rambaut Exp $
 */
public class BranchRatePlotter {

    public static void main(String[] args) throws java.io.IOException, Importer.ImportException {

        String controlFile = args[0];

        //String treeFile1 = args[0];
        //String treeFile2 = args[1];

        String targetTreeFile = args[1];

        int burnin = 0;
        if (args.length > 2) {
            burnin = Integer.parseInt(args[2]);
        }
        System.out.println("Ignoring first " + burnin + " trees as burnin.");

        BufferedReader readerTarget = new BufferedReader(new FileReader(targetTreeFile));
        String lineTarget = readerTarget.readLine();
        readerTarget.close();

        TreeImporter targetImporter;
        if (lineTarget.toUpperCase().startsWith("#NEXUS")) {
            targetImporter = new NexusImporter(new FileReader(targetTreeFile));
        } else {
            targetImporter = new NewickImporter(new FileReader(targetTreeFile));
        }
        MutableTree targetTree = new FlexibleTree(targetImporter.importNextTree());
        targetTree = TreeUtils.rotateTreeByComparator(targetTree, TreeUtils.createNodeDensityComparator(targetTree));

        BufferedReader reader = new BufferedReader(new FileReader(controlFile));
        String line = reader.readLine();

        int totalTrees = 0;
        int totalTreesUsed = 0;
        while (line != null) {

            StringTokenizer tokens = new StringTokenizer(line);

            NexusImporter importer1 = new NexusImporter(new FileReader(tokens.nextToken()));
            NexusImporter importer2 = new NexusImporter(new FileReader(tokens.nextToken()));

            int fileTotalTrees = 0;
            while (importer1.hasTree()) {
                Tree timeTree = importer1.importNextTree();
                Tree mutationTree = importer2.importNextTree();

                if (fileTotalTrees >= burnin) {
                    annotateRates(targetTree, targetTree.getRoot(), timeTree, mutationTree);

                    totalTreesUsed += 1;
                }
                totalTrees += 1;
                fileTotalTrees += 1;
            }
            line = reader.readLine();
        }

        System.out.println("Total trees read: " + totalTrees);
        System.out.println("Total trees summarized: " + totalTreesUsed);

        // collect all rates
        double mutations = 0.0;
        double time = 0.0;

        double[] rates = new double[targetTree.getNodeCount()-1];
        int index = 0;
        for (int i = 0; i < targetTree.getNodeCount(); i++) {
            NodeRef node = targetTree.getNode(i);

            if (!targetTree.isRoot(node)) {
	            Integer count = ((Integer)targetTree.getNodeAttribute(node,"count"));
	            if (count == null) {
		            throw new RuntimeException("Count missing from node in target tree");
	            }
                if (!targetTree.isExternal(node)) {
                    double prob = (double) (int) count /(double)(totalTreesUsed);
                    if (prob >= 0.5) {
                        String label = ""+(Math.round(prob*100)/100.0);
                        targetTree.setNodeAttribute(node, "label", label);
                    }
                }

                Number totalMutations = (Number)targetTree.getNodeAttribute(node,"totalMutations");
                Number totalTime = (Number)targetTree.getNodeAttribute(node,"totalTime");
                mutations += totalMutations.doubleValue();
                time += totalTime.doubleValue();

                rates[index] = totalMutations.doubleValue()/totalTime.doubleValue();
                System.out.println(totalMutations.doubleValue() + " / " + totalTime.doubleValue() + " = " + rates[index]);
                targetTree.setNodeRate(node, rates[index]);
                index += 1;
            }
        }

        double minRate = DiscreteStatistics.min(rates);
        double maxRate = DiscreteStatistics.max(rates);
        double medianRate = DiscreteStatistics.median(rates);
        //double topThird = DiscreteStatistics.quantile(2.0/3.0,rates);
        //double bottomThird = DiscreteStatistics.quantile(1.0/3.0,rates);
        //double unweightedMeanRate = DiscreteStatistics.mean(rates);
        double meanRate = mutations/time;
        System.out.println(minRate + "\t" + maxRate + "\t" + medianRate + "\t" + meanRate);

        for (int i = 0; i < targetTree.getNodeCount(); i++) {
            NodeRef node = targetTree.getNode(i);

            if (!targetTree.isRoot(node)) {
                double rate = targetTree.getNodeRate(node);

                //double branchTime = ((Number)targetTree.getNodeAttribute(node, "totalTime")).doubleValue();
                //double branchMutations = ((Number)targetTree.getNodeAttribute(node, "totalMutations")).doubleValue();

                float relativeRate = (float)(rate /maxRate);
                float radius = (float)Math.sqrt(relativeRate*36.0);

                /*
                float relativeRateZero = (float)((rate - minRate) / (maxRate - minRate));

                float red = 0.0f;
                float green = 0.0f;
                float blue = 0.0f;
                if (relativeRateZero < 0.5f) {
                    blue =  1.0f - (relativeRateZero * 2.0f);
                    green = 1.0f - blue;
                } else {
                    red = (relativeRateZero - 0.5f) * 2.0f;
                    green = 1.0f - red;
                }
                */
                //System.out.println(red + " " + green + " " + blue);

                //float lineThickness = relativeRate*6.0f;

                if (rate > meanRate) {
                    targetTree.setNodeAttribute(node, "color", new Color(1.0f, 0.5f, 0.5f));
                } else {
                    targetTree.setNodeAttribute(node, "color", new Color(0.5f, 0.5f, 1.0f));
                }
                //targetTree.setNodeAttribute(node, "color", new Color(red, green, blue));
                targetTree.setNodeAttribute(node, "line", new BasicStroke(1.0f));
                targetTree.setNodeAttribute(node, "shape",
                        new java.awt.geom.Ellipse2D.Double(0,0,radius*2.0,radius*2.0));

            }

            java.util.List heightList = (java.util.List)targetTree.getNodeAttribute(node, "heightList");
            if (heightList != null) {
                double[] heights = new double[heightList.size()];
                for (int j = 0; j < heights.length; j++) {
                    heights[j] = (Double) heightList.get(j);
                }
                targetTree.setNodeHeight(node, DiscreteStatistics.mean(heights));
                //if (heights.length >= (totalTreesUsed/2)) {
                    targetTree.setNodeAttribute(node, "nodeHeight.mean", DiscreteStatistics.mean(heights));
                    targetTree.setNodeAttribute(node, "nodeHeight.hpdUpper", DiscreteStatistics.quantile(0.975, heights));
                    targetTree.setNodeAttribute(node, "nodeHeight.hpdLower", DiscreteStatistics.quantile(0.025, heights));
                    //targetTree.setNodeAttribute(node, "nodeHeight.max", new Double(DiscreteStatistics.max(heights)));
                    //targetTree.setNodeAttribute(node, "nodeHeight.min", new Double(DiscreteStatistics.min(heights)));
                //}
            }
        }

        StringBuffer buffer = new StringBuffer();
        writeTree(targetTree, targetTree.getRoot(), buffer, true, false);
        buffer.append(";\n");
        writeTree(targetTree, targetTree.getRoot(), buffer, false, true);
        buffer.append(";\n");
        System.out.println(buffer.toString());

        SquareTreePainter treePainter = new SquareTreePainter();
        treePainter.setColorAttribute("color");
        treePainter.setLineAttribute("line");
//        treePainter.setShapeAttribute("shape");
//        treePainter.setLabelAttribute("label");
        JTreeDisplay treeDisplay = new JTreeDisplay(treePainter,targetTree);

        JTreePanel treePanel = new JTreePanel(treeDisplay);

        JFrame frame = new JFrame();
        frame.setSize(800,600);
        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(treePanel);
        frame.setVisible(true);

        PrinterJob printJob = PrinterJob.getPrinterJob();
        printJob.setPrintable(treeDisplay);
        if (printJob.printDialog()){
            try{
                  printJob.print();
            }
            catch(Exception ex){
                  throw new RuntimeException(ex);
            }
        }
    }

    private static void writeTree(Tree tree, NodeRef node, StringBuffer buffer, boolean rates, boolean labels) {

        if (tree.isExternal(node)) {
		    buffer.append(tree.getTaxonId(node.getNumber()));

        } else {
            buffer.append("(");
            writeTree(tree, tree.getChild(node, 0), buffer, rates, labels);
            for (int i = 1; i < tree.getChildCount(node); i++) {
                buffer.append(",");
                writeTree(tree, tree.getChild(node, i), buffer, rates, labels);
            }
            buffer.append(")");
        }

        NodeRef parent = tree.getParent(node);
        if (parent != null) {
            double totalMutations = (Double) tree.getNodeAttribute(node, "totalMutations");
            double totalTime = (Double) tree.getNodeAttribute(node, "totalTime");
            double rate = totalMutations/totalTime;

            int count = (Integer) tree.getNodeAttribute(node, "count");
            if (rates) {
                buffer.append(":").append(String.valueOf(rate));
            } else {
                buffer.append(":").append(String.valueOf(count));
            }
        }
    }

    private static void annotateRates(
            MutableTree targetTree, NodeRef node, Tree timeTree, Tree mutationTree) {

        Set<String> leafSet = TreeUtils.getDescendantLeaves(targetTree, node);
        if (TreeUtils.isMonophyletic(timeTree, leafSet)) {
            NodeRef timeNode = TreeUtils.getCommonAncestorNode(timeTree, leafSet);
            NodeRef mutationNode = TreeUtils.getCommonAncestorNode(mutationTree, leafSet);

            double height = timeTree.getNodeHeight(timeNode);

            if (!targetTree.isRoot(node)) {
                double time = timeTree.getNodeHeight(timeTree.getParent(timeNode)) - height;
                double mutations = mutationTree.getNodeHeight(mutationTree.getParent(mutationNode)) - mutationTree.getNodeHeight(mutationNode);

                //double rate = mutations/time;
                Number totalMutations = (Number)targetTree.getNodeAttribute(node, "totalMutations");
                Number totalTime = (Number)targetTree.getNodeAttribute(node, "totalTime");
                if (totalMutations == null) {
                    targetTree.setNodeAttribute(node, "totalMutations", mutations);
                    targetTree.setNodeAttribute(node, "totalTime", time);
                    targetTree.setNodeAttribute(node, "count", 1);
                } else {
                    Integer count = (Integer)targetTree.getNodeAttribute(node, "count");
                    targetTree.setNodeAttribute(node, "totalMutations", totalMutations.doubleValue() + mutations);
                    targetTree.setNodeAttribute(node, "totalTime", totalTime.doubleValue() + time);
                    targetTree.setNodeAttribute(node, "count", count + 1);
                }
            }
            if (!targetTree.isExternal(node)) {
                java.util.List<Double> list = (java.util.List<Double>)targetTree.getNodeAttribute(node, "heightList");
                if (list == null) {
                    list = new ArrayList<Double>() ;
                    targetTree.setNodeAttribute(node, "heightList", list);
                }
                list.add(height);
            }
        }


        for (int i = 0; i < targetTree.getChildCount(node); i++) {
            annotateRates(targetTree, targetTree.getChild(node, i), timeTree, mutationTree);
        }
    }
}
