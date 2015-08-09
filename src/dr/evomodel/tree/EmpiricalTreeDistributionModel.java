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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.math.MathUtils;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.model.Statistic;
import dr.inference.model.StatisticList;
import dr.inference.model.AbstractModel;

import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class EmpiricalTreeDistributionModel extends TreeModel {

    public EmpiricalTreeDistributionModel(final Tree[] trees, int startingTree) {
        super(EMPIRICAL_TREE_DISTRIBUTION_MODEL);

        this.trees = trees;
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
        storedTreeIndex = currentTreeIndex;
    }

    protected void restoreState() {
        currentTreeIndex = storedTreeIndex;
    }

    protected void acceptState() {
    }

    public void drawTreeIndex() {
        drawTreeIndex(-1);
    }

    private void drawTreeIndex(int treeNumber) {
//        System.err.print("Drawing new tree, (old tree = " + currentTreeIndex);

        if (treeNumber == -1) {
            currentTreeIndex = MathUtils.nextInt(trees.length);
        } else {
            currentTreeIndex = treeNumber;
        }

        // Force computation of node heights now rather than later in the evaluation
        // where multithreading may get conflicts.
        trees[currentTreeIndex].getNodeHeight(trees[currentTreeIndex].getRoot());

//        System.err.println(") new tree = " + currentTreeIndex);

        fireModelChanged(new TreeModel.TreeChangedEvent());
    }

    public NodeRef getRoot() {
        return trees[currentTreeIndex].getRoot();
    }

    public int getNodeCount() {
        return trees[currentTreeIndex].getNodeCount();
    }

    public NodeRef getNode(final int i) {
        return trees[currentTreeIndex].getNode(i);
    }

    public NodeRef getInternalNode(final int i) {
        return trees[currentTreeIndex].getInternalNode(i);
    }

    public NodeRef getExternalNode(final int i) {
        return trees[currentTreeIndex].getExternalNode(i);
    }

    public int getExternalNodeCount() {
        return trees[currentTreeIndex].getExternalNodeCount();
    }

    public int getInternalNodeCount() {
        return trees[currentTreeIndex].getInternalNodeCount();
    }

    public Taxon getNodeTaxon(final NodeRef node) {
        return trees[currentTreeIndex].getNodeTaxon(node);
    }

    public boolean hasNodeHeights() {
        return trees[currentTreeIndex].hasNodeHeights();
    }

    public double getNodeHeight(final NodeRef node) {
        return trees[currentTreeIndex].getNodeHeight(node);
    }

    public boolean hasBranchLengths() {
        return trees[currentTreeIndex].hasBranchLengths();
    }

    public double getBranchLength(final NodeRef node) {
        return trees[currentTreeIndex].getBranchLength(node);
    }

    public double getNodeRate(final NodeRef node) {
        return trees[currentTreeIndex].getNodeRate(node);
    }

    public Object getNodeAttribute(final NodeRef node, final String name) {
        return trees[currentTreeIndex].getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(final NodeRef node) {
        return trees[currentTreeIndex].getNodeAttributeNames(node);
    }

    public boolean isExternal(final NodeRef node) {
        return trees[currentTreeIndex].isExternal(node);
    }

    public boolean isRoot(final NodeRef node) {
        return trees[currentTreeIndex].isRoot(node);
    }

    public int getChildCount(final NodeRef node) {
        return trees[currentTreeIndex].getChildCount(node);
    }

    public NodeRef getChild(final NodeRef node, final int j) {
        return trees[currentTreeIndex].getChild(node, j);
    }

    public NodeRef getParent(final NodeRef node) {
        return trees[currentTreeIndex].getParent(node);
    }

    public Tree getCopy() {
        return trees[currentTreeIndex].getCopy();
    }

    public int getTaxonCount() {
        return trees[currentTreeIndex].getTaxonCount();
    }

    public Taxon getTaxon(final int taxonIndex) {
        return trees[currentTreeIndex].getTaxon(taxonIndex);
    }

    public String getTaxonId(final int taxonIndex) {
        return trees[currentTreeIndex].getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(final String id) {
        return trees[currentTreeIndex].getTaxonIndex(id);
    }

    public int getTaxonIndex(final Taxon taxon) {
        return trees[currentTreeIndex].getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        return trees[currentTreeIndex].asList();
    }

    public Object getTaxonAttribute(final int taxonIndex, final String name) {
        return trees[currentTreeIndex].getTaxonAttribute(taxonIndex, name);
    }

    public Iterator<Taxon> iterator() {
        return trees[currentTreeIndex].iterator();
    }

    public Type getUnits() {
        return trees[currentTreeIndex].getUnits();
    }

    public void setUnits(final Type units) {
        trees[currentTreeIndex].setUnits(units);
    }

    public void setAttribute(final String name, final Object value) {
        trees[currentTreeIndex].setAttribute(name, value);
    }

    public Object getAttribute(final String name) {
        return trees[currentTreeIndex].getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return trees[currentTreeIndex].getAttributeNames();
    }

    public static final String EMPIRICAL_TREE_DISTRIBUTION_MODEL = "empiricalTreeDistributionModel";

    private final Tree[] trees;
    private int currentTreeIndex;
    private int storedTreeIndex;
}
