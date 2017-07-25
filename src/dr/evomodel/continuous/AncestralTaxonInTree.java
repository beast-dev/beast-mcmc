/*
 * RestrictedPartials.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class AncestralTaxonInTree extends AbstractModel {

    public AncestralTaxonInTree(Taxon ancestor,
                                MutableTreeModel treeModel,
                                TaxonList taxonList,
                                Parameter priorSampleSize) throws TreeUtils.MissingTaxonException {
        this(ancestor, treeModel, taxonList, priorSampleSize, null, -1);
    }

    public AncestralTaxonInTree(Taxon ancestor,
                                MutableTreeModel treeModel,
                                TaxonList descendents,
                                Parameter priorSampleSize,
                                NodeRef node, int index) throws TreeUtils.MissingTaxonException {

        super(ancestor.getId());

        this.ancestor = ancestor;
        this.treeModel = treeModel;
        this.descendents = descendents;
        this.pseudoBranchLength = priorSampleSize;
        this.index = index;
        this.node = node;

        this.tips = TreeUtils.getTipsForTaxa(treeModel, descendents);
        this.tipBitSet = TreeUtils.getTipsBitSetForTaxa(treeModel, descendents);

        addVariable(priorSampleSize);
    }

    // Public API

    final public double getPseudoBranchLength() {
        return pseudoBranchLength.getParameterValue(0);
    }

    final MutableTreeModel getTreeModel() { return treeModel; }

//    final double[] getPartials() { return meanParameter.getParameterValues(); }

//    final double getPartial(int i) { return meanParameter.getParameterValue(i); }

//    final double getPriorSampleSize() { return pseudoBranchLength.getParameterValue(0); }

    final double[] getRestrictedPartials() {
        assert(false);
        return null;
    }

//    final PrecisionType getPrecisionType() {
//        assert(false);
//        return PrecisionType.SCALAR;
//    }

    final public int getIndex() { return index; }

    final public void setIndex(int index) { this.index = index; }

    final public NodeRef getNode() { return node; }

    final public void setNode(NodeRef node) { this.node = node; }

    // AbstractModel implementation

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged(variable, index);
    }

    final private Taxon ancestor;
    final private MutableTreeModel treeModel;
    final private TaxonList descendents;

    final private Set<Integer> tips;
    final private BitSet tipBitSet;

    final private Parameter pseudoBranchLength;

    private int index;
    private NodeRef node;

    public TaxonList getTaxonList() {
        return descendents;
    }

    public BitSet getTipBitSet() {
        return tipBitSet;
    }

    public Taxon getTaxon() {
        return ancestor;
    }
}
