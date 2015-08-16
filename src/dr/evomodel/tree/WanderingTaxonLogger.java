/*
 * WanderingTaxonLogger.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class WanderingTaxonLogger implements TreeTraitProvider {

    public static final String WANDERER = "wanderingTaxonLogger";
    public static final String RELATIVE = "relative";

    public WanderingTaxonLogger(String name, Taxon taxon, Relative relative) {
        if (name == null) {
            this.name = RELATIVE;
        } else {
            this.name = name;
        }
        this.taxon = taxon;
        this.relative = relative;
    }

    TreeTrait relativeTrait = new TreeTrait.I() {
        public String getTraitName() {
            return name;
        }

        public Intent getIntent() {
            return Intent.NODE;
        }

        public Integer getTrait(Tree tree, NodeRef node) {
            int rtnValue = 0;
            if (relative == Relative.PARENT) {
                if (isAnyChildEqualToTaxon(tree, node, taxon, null)) {
                    rtnValue = 1;
                }
            } else if (relative == Relative.SISTER && !tree.isRoot(node)) {
                if (isAnyChildEqualToTaxon(tree, tree.getParent(node), taxon, node)) {
                    rtnValue = 1;
                }
            }

            return rtnValue;
        }
    };

    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[] { relativeTrait };
    }

    public TreeTrait getTreeTrait(String key) {
        // ignore the key - it must be the one they wanted, no?
        return relativeTrait;
    }

    private boolean isAnyChildEqualToTaxon(Tree tree, NodeRef node, Taxon taxon, NodeRef exclude) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (child != exclude && tree.isExternal(child)) {
                String taxonString = tree.getNodeTaxon(child).getId();
                if (taxonString.equals(taxon.getId())) {
                    return true;
                }
            }
        }
        return false;
    }

    public enum Relative {
        PARENT,
        SISTER
    }

    private String name;
    private Taxon taxon;
    private Relative relative;

}
