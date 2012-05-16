/*
 * TreeTransform.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public abstract class TreeTransform extends AbstractModel implements TreeTraitProvider, Citable {

    public static final String TREE_TRANSFORM_PREFIX = "treeTransform";

    public TreeTransform(String name) {
        super(name);
        setupTraits();
    }

    private void setupTraits() {
        TreeTrait baseTrait = new TreeTrait.D() {

            public String getTraitName() {
                return TREE_TRANSFORM_PREFIX;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public Double getTrait(Tree tree, NodeRef node) {
                return getScaleForNode(tree, node);
            }

            public boolean getLoggable() {
                return true;
            }
        };
        treeTraits.addTrait(baseTrait);
    }

    public abstract double transform(TransformedTreeModel tree, NodeRef node, double originalHeight);

    protected abstract double getScaleForNode(Tree tree, NodeRef node);

    public abstract String getInfo();

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    private final Helper treeTraits = new Helper();

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("P", "Lemey"),
                                new Author("MA", "Suchard"),
                        },
                        Citation.Status.IN_PREPARATION
                ));
        return citations;
    }
}

