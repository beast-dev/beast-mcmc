/*
 * TreeTransform.java
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
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

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

    // TODO originalHeight is not necessary to pass anymore; remove
    public abstract double transform(Tree tree, NodeRef node, double originalHeight);

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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Branch-specific phenotypic mixture model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.VRANCKEN_2015_SIMULTANEOUSLY);
    }


}

