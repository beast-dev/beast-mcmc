/*
 * TransformedTreeModel.java
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

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;
import java.util.logging.Logger;

/**
 * A transformation of a tree model
 *
 * @author Marc Suchard
 */
public class TransformedTreeModel extends AbstractModel implements MutableTreeModel, TransformableTree, Citable {

    public TransformedTreeModel(String id, MutableTreeModel tree, TreeTransform treeTransform) {
        super(id);
        this.treeModel = tree;
        this.treeTransform = treeTransform;
        addModel(tree);
        addModel(treeTransform);

        Logger log = Logger.getLogger("dr.evomodel.tree");
        log.info("Creating a transform tree.");
    }

    public String toString() {
        return TreeUtils.newick(this);
    }

    public double getNodeHeight(NodeRef node) {
        return treeTransform.transform(treeModel, node, treeModel.getNodeHeight(node));
    }

    public double getBranchLength(NodeRef node) {
        NodeRef parent = treeModel.getParent(node);
        if (parent == null) {
            return 0.0;
        }
//        System.err.println("p: " + this.getNodeHeight(parent));
//        System.err.println("c: " + this.getNodeHeight(node));
        return this.getNodeHeight(parent) - this.getNodeHeight(node);
    }

//    public double getOriginalNodeHeight(NodeRef node) {
//        return treeModel.getNodeHeight(node);
//    }

    // TODO 1. Reparameterize via parentNodeHeight
    // TODO 2. Deal with branchRateModel issues and no values from external nodes

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeTransform) {
            fireModelChanged(new TreeChangedEvent.WholeTree()); // All internal node heights have changed!
        } else if (model == treeModel) {
            fireModelChanged(object, index);
        } else {
            throw new IllegalArgumentException("Illegal model");
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // Do nothing; no variables
    }

    // Delegate the rest to treeModel

    public NodeRef getRoot() {
        return treeModel.getRoot();
    }

    public int getNodeCount() {
        return treeModel.getNodeCount();
    }

    public NodeRef getNode(int i) {
        return treeModel.getNode(i);
    }

    public NodeRef getInternalNode(int i) {
        return treeModel.getInternalNode(i);
    }

    public NodeRef getExternalNode(int i) {
        return treeModel.getExternalNode(i);
    }

    public int getExternalNodeCount() {
        return treeModel.getExternalNodeCount();
    }

    public int getInternalNodeCount() {
        return treeModel.getInternalNodeCount();
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return treeModel.getNodeTaxon(node);
    }

    public boolean hasNodeHeights() {
        return treeModel.hasNodeHeights();
    }

    public boolean hasBranchLengths() {
        return treeModel.hasBranchLengths();
    }

    public double getNodeRate(NodeRef node) {
        return treeModel.getNodeRate(node);
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        return treeModel.getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        return treeModel.getNodeAttributeNames(node);
    }

    public boolean isExternal(NodeRef node) {
        return treeModel.isExternal(node);
    }

    public boolean isRoot(NodeRef node) {
        return treeModel.isRoot(node);
    }

    public int getChildCount(NodeRef node) {
        return treeModel.getChildCount(node);
    }

    public NodeRef getChild(NodeRef node, int j) {
        return treeModel.getChild(node, j);
    }

    public NodeRef getParent(NodeRef node) {
        return treeModel.getParent(node);
    }

    public Tree getCopy() {
        return treeModel.getCopy();
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        // Do nothing
    }

    protected void acceptState() {
        // Do nothing
    }

    private final TreeTransform treeTransform;
    private final MutableTreeModel treeModel;

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        return treeModel.getMultivariateNodeTrait(node, name);
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        treeModel.setMultivariateTrait(n, name, value);
    }

//    public Parameter getRootHeightParameter() {
//        return treeModel.getRootHeightParameter();
//    }

    public boolean beginTreeEdit() {
        return treeModel.beginTreeEdit();
    }

    public void endTreeEdit() {
        treeModel.endTreeEdit();
    }

    public void addChild(NodeRef parent, NodeRef child) {
        treeModel.addChild(parent, child);
    }

    public void removeChild(NodeRef parent, NodeRef child) {
        treeModel.removeChild(parent, child);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        treeModel.replaceChild(node, child, newChild);
    }

    public void setRoot(NodeRef root) {
        treeModel.setRoot(root);
    }

    public void setNodeHeight(NodeRef node, double height) {
        treeModel.setNodeHeight(node, height);
    }

    public void setNodeRate(NodeRef node, double height) {
       treeModel.setNodeRate(node, height);
    }

    public void setBranchLength(NodeRef node, double length) {
        treeModel.setBranchLength(node, length);
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        treeModel.setNodeAttribute(node, name, value);
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        treeModel.addMutableTreeListener(listener);
    }

    public void setAttribute(String name, Object value) {
        treeModel.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        return treeModel.getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return treeModel.getAttributeNames();
    }

    public int addTaxon(Taxon taxon) {
        return treeModel.addTaxon(taxon);
    }

    public boolean removeTaxon(Taxon taxon) {
        return treeModel.removeTaxon(taxon);
    }

    public void setTaxonId(int taxonIndex, String id) {
        treeModel.setTaxonId(taxonIndex, id);
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        treeModel.setTaxonAttribute(taxonIndex, name, value);
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        treeModel.addMutableTaxonListListener(listener);
    }

    public int getTaxonCount() {
        return treeModel.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex) {
        return treeModel.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex) {
        return treeModel.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id) {
        return treeModel.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon) {
        return treeModel.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        return treeModel.asList();
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        return treeModel.getTaxonAttribute(taxonIndex, name);
    }

    public Iterator<Taxon> iterator() {
        return treeModel.iterator();
    }

    public Type getUnits() {
        return treeModel.getUnits();
    }

    public void setUnits(Type units) {
        treeModel.setUnits(units);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian estimation of Pagel's lambda";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.VRANCKEN_2015_SIMULTANEOUSLY);
    }

    @Override
    public NodeRef getOriginalNode(NodeRef transformedNode) {
        return transformedNode;
    }

    @Override
    public NodeRef getTransformedNode(NodeRef originalNode) {
        return originalNode;
    }

    @Override
    public Tree getOriginalTree() {
        return treeModel;
    }
}
