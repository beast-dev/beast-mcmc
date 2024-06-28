/*
 * TransformedTreeTraitProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.tree;

import dr.util.Transform;

/**
 * @author Alexander Fisher
 * @author Marc A. Suchard
 */
public class TransformedTreeTraitProvider implements TreeTraitProvider {

    public TransformedTreeTraitProvider(TreeTraitProvider treeTraitProvider,
                                        Transform transform) {

        for (TreeTrait trait : treeTraitProvider.getTreeTraits()) {
            if (trait.getTraitClass() == Double.class) {
                transformedTreeTraits.addTrait(createD(trait, transform));
            } else if (trait.getTraitClass() == double[].class) {
                transformedTreeTraits.addTrait(createDA(trait, transform));
            } else {
                throw new IllegalArgumentException("Not transformable trait: " + trait.getTraitName());
            }
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return transformedTreeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return transformedTreeTraits.getTreeTrait(key);
    }

    private TreeTrait.D createD(final TreeTrait originalTrait, final Transform transform) {

        return new TreeTrait.D() {
            @Override
            public String getTraitName() {
                return transform.getTransformName() + "." + originalTrait.getTraitName();
            }

            @Override
            public Intent getIntent() {
                return originalTrait.getIntent();
            }

            @Override
            public Double getTrait(Tree tree, NodeRef node) {
                return transform.transform((Double) originalTrait.getTrait(tree, node));
            }

            @Override
            public boolean getLoggable() {
                return true;
            }
        };
    }

    private TreeTrait.DA createDA(final TreeTrait originalTrait, final Transform transform) {

        return new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return transform.getTransformName() + "." + originalTrait.getTraitName();
            }

            @Override
            public Intent getIntent() {
                return originalTrait.getIntent();
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                double[] raw = (double[]) originalTrait.getTrait(tree, node);
                return transform.transform(raw, 0, raw.length);
            }

            @Override
            public boolean getLoggable() {
                return true;
            }
        };
    }

    private Helper transformedTreeTraits = new Helper();
}
