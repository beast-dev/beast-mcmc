/*
 * CladeRelationshipStatistic.java
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
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.inference.model.BooleanStatistic;
import dr.math.UnivariateMinimum;

import java.util.Collections;
import java.util.Set;

/**
 * Assesses relationship between pairs of clades.
 *
 * @author Andy Magee
 */
public class CladeRelationshipStatistic extends TreeStatistic implements BooleanStatistic {

    public CladeRelationshipStatistic(String name,
                                      Tree tree,
                                      TaxonList taxaA,
                                      TaxonList taxaB,
                                      RelationshipType type) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;
        this.leafSetA = TreeUtils.getLeavesForTaxa(tree, taxaA);
        this.leafSetB = TreeUtils.getLeavesForTaxa(tree, taxaB);

        this.type = type;
    }

    public enum RelationshipType {
        SISTER("sister") {
            public boolean extractResultForType(boolean[] results) {
                return results[0];
            }
        },
        A_IN_B("aInB") {
            public boolean extractResultForType(boolean[] results) {
                return results[1];
            }
        };

        RelationshipType(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public abstract boolean extractResultForType(boolean[] results);

        private String name;
    }


    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    public boolean getBoolean(int dim) {

        boolean[] results = getRelationship();

        return type.extractResultForType(results);

   }

    public boolean[] getRelationship() {
        boolean[] results = new boolean[2];

        NodeRef mrcaA = TreeUtils.getCommonAncestorNode(tree, leafSetA);
        NodeRef mrcaB = TreeUtils.getCommonAncestorNode(tree, leafSetB);
        NodeRef mrcaTotal = TreeUtils.getCommonAncestorNode(tree, mrcaA, mrcaB);

        // Sister
        if ( tree.getParent(mrcaA) == tree.getParent(mrcaB)) {
            results[0] = true;
        } else if (mrcaTotal == mrcaB) { // A in B
            results[1] = true;
        }

        return results;
    }

    public static RelationshipType factory(String match) {
        for (CladeRelationshipStatistic.RelationshipType type : CladeRelationshipStatistic.RelationshipType.values()) {
            if (match.equalsIgnoreCase(type.name)) {
                return type;
            }
        }
        return null;
    }


    private Tree tree;
    private final Set<String> leafSetA;
    private final Set<String> leafSetB;
    private final RelationshipType type;
}
