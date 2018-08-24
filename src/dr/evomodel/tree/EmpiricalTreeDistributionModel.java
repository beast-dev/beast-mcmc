/*
 * EmpiricalTreeDistributionModel.java
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

package dr.evomodel.tree;

import dr.evolution.io.Importer;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;
import dr.inference.model.Statistic;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class EmpiricalTreeDistributionModel extends TreeModel {

    /**
     * This constructor takes an array of trees and jumps randomly amongst them.
     * @param trees
     * @param startingTree
     */
    public EmpiricalTreeDistributionModel(final Tree[] trees, int startingTree) {
        this(trees, null, startingTree);
    }

    /**
     * This constructor takes a TreeImporter and reads the trees one by one to the end (and
     * then starts throwing exceptions).
     * @param importer
     * @param startingTree
     */
    public EmpiricalTreeDistributionModel(final TreeImporter importer, int startingTree) {
        this(null, importer, startingTree);
    }

    private EmpiricalTreeDistributionModel(final Tree[] trees, final TreeImporter importer, int startingTree) {
        super(EMPIRICAL_TREE_DISTRIBUTION_MODEL);

        this.trees = trees;
        this.importer = importer;
        drawTreeIndex(startingTree);

        addStatistic(new Statistic.Abstract("Current Tree")  {

            public int getDimension() {
                return 1;
            }

            public double getStatisticValue(int dim) {
                return currentTreeIndex;
            }
        });
    }

    protected void storeState() {
        storedCurrentTree = currentTree;
    }

    protected void restoreState() {
        currentTree = storedCurrentTree;
    }

    protected void acceptState() {
    }

    public void drawTreeIndex() {
        drawTreeIndex(-1);
    }

    private void drawTreeIndex(int treeNumber) {
//        System.err.print("Drawing new tree, (old tree = " + currentTreeIndex);

        if (importer != null) {
            try {
                if (importer.hasTree() == false) {
                    throw new RuntimeException("EmpiricalTreeDistributionModel has run out of trees");
                }
                currentTree = importer.importNextTree();
            } catch (IOException e) {
                throw new RuntimeException("EmpiricalTreeDistributionModel unable to load next tree");
            } catch (Importer.ImportException e) {
                throw new RuntimeException("EmpiricalTreeDistributionModel unable to load next tree");
            }
            currentTreeIndex += 1;
        } else {
            if (treeNumber == -1) {
                currentTreeIndex = MathUtils.nextInt(trees.length);
                currentTree = trees[currentTreeIndex];
            } else {
                currentTreeIndex = treeNumber;
                currentTree = trees[currentTreeIndex];
            }
        }

        // Force computation of node heights now rather than later in the evaluation
        // where multithreading may get conflicts.
        currentTree.getNodeHeight(currentTree.getRoot());

//        System.err.println(") new tree = " + currentTreeIndex);

        fireModelChanged(new TreeChangedEvent());
    }

    public NodeRef getRoot() {
        return currentTree.getRoot();
    }

    public int getNodeCount() {
        return currentTree.getNodeCount();
    }

    public NodeRef getNode(final int i) {
        return currentTree.getNode(i);
    }

    public NodeRef getInternalNode(final int i) {
        return currentTree.getInternalNode(i);
    }

    public NodeRef getExternalNode(final int i) {
        return currentTree.getExternalNode(i);
    }

    public int getExternalNodeCount() {
        return currentTree.getExternalNodeCount();
    }

    public int getInternalNodeCount() {
        return currentTree.getInternalNodeCount();
    }

    public Taxon getNodeTaxon(final NodeRef node) {
        return trees[currentTreeIndex].getNodeTaxon(node);
    }

    public boolean hasNodeHeights() {
        return currentTree.hasNodeHeights();
    }

    public double getNodeHeight(final NodeRef node) {
        return currentTree.getNodeHeight(node);
    }

    public boolean hasBranchLengths() {
        return currentTree.hasBranchLengths();
    }

    public double getBranchLength(final NodeRef node) {
        return currentTree.getBranchLength(node);
    }

    public double getNodeRate(final NodeRef node) {
        return currentTree.getNodeRate(node);
    }

    public Object getNodeAttribute(final NodeRef node, final String name) {
        return currentTree.getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(final NodeRef node) {
        return currentTree.getNodeAttributeNames(node);
    }

    public boolean isExternal(final NodeRef node) {
        return currentTree.isExternal(node);
    }

    public boolean isRoot(final NodeRef node) {
        return currentTree.isRoot(node);
    }

    public int getChildCount(final NodeRef node) {
        return currentTree.getChildCount(node);
    }

    public NodeRef getChild(final NodeRef node, final int j) {
        return currentTree.getChild(node, j);
    }

    public NodeRef getParent(final NodeRef node) {
        return currentTree.getParent(node);
    }

    public Tree getCopy() {
        return currentTree.getCopy();
    }

    public int getTaxonCount() {
        return currentTree.getTaxonCount();
    }

    public Taxon getTaxon(final int taxonIndex) {
        return currentTree.getTaxon(taxonIndex);
    }

    public String getTaxonId(final int taxonIndex) {
        return currentTree.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(final String id) {
        return currentTree.getTaxonIndex(id);
    }

    public int getTaxonIndex(final Taxon taxon) {
        return currentTree.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        return currentTree.asList();
    }

    public Object getTaxonAttribute(final int taxonIndex, final String name) {
        return currentTree.getTaxonAttribute(taxonIndex, name);
    }

    public Iterator<Taxon> iterator() {
        return currentTree.iterator();
    }

    public Type getUnits() {
        return currentTree.getUnits();
    }

    public void setUnits(final Type units) {
        currentTree.setUnits(units);
    }

    public void setAttribute(final String name, final Object value) {
        currentTree.setAttribute(name, value);
    }

    public Object getAttribute(final String name) {
        return currentTree.getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return currentTree.getAttributeNames();
    }

    public static final String EMPIRICAL_TREE_DISTRIBUTION_MODEL = "empiricalTreeDistributionModel";

    private final Tree[] trees;
    private final TreeImporter importer;
    private Tree currentTree;
    private Tree storedCurrentTree;

    private int currentTreeIndex;
}
