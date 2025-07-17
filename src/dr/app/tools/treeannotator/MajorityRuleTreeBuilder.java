/*
 * HIPSTRTreeBuilder.java
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

package dr.app.tools.treeannotator;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.*;

public class MajorityRuleTreeBuilder {

    public MutableTree getMajorityRuleConsensusTree(CladeSystem cladeSystem, TaxonList taxonList) {
        BiClade rootClade = (BiClade)cladeSystem.getRootClade();

        for (BiClade mrClade : cladeSystem.getTipClades()) {
            BiClade mrParent = findMajorityRuleParent(mrClade);
            mrParent.addChild(mrClade);
        }
        Set<BiClade> majorityRuleClades = cladeSystem.getTopClades(0.5);
        for (BiClade mrClade : majorityRuleClades) {
            if (mrClade != rootClade) {
                BiClade mrParent = findMajorityRuleParent(mrClade);
                mrParent.addChild(mrClade);
                assert true;
            }
        }

        // create a map so that tip numbers are in the same order as the taxon list
        Map<Taxon, Integer> taxonNumberMap = new HashMap<>();
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonNumberMap.put(taxonList.getTaxon(i), i);
        }

        FlexibleNode root = buildMajorityRuleConsensusTree(rootClade);
        FlexibleTree tree = new FlexibleTree(root, taxonNumberMap);
        tree.setLengthsKnown(false);
        return tree;
    }

//    private static BiClade findMajorityRuleParent(BiClade clade) {
//        if (clade.getMajorityRuleParent() != null) {
//            return clade.getMajorityRuleParent();
//        }
//
//        assert clade.getParentClades() != null;
//
//        BiClade mrParent = null;
//
//        // find the smallest parent clade that is majority rule
//        for (BiClade parent : clade.getParentClades()) {
//            if (parent.getCredibility() >= 0.5) {
//                assert mrParent == null || parent.getSize() != mrParent.getSize();
//                if (mrParent == null || parent.getSize() < mrParent.getSize()) {
//                    mrParent = parent;
//                }
//            }
//        }
//
//        if (mrParent == null) {
//            // no MR direct parents so crawl upwards
//            for (BiClade parent : clade.getParentClades()) {
//                BiClade p = parent;
//                while (p.getCredibility() < 0.5 && p.getParentClades() != null) {
//                    p = findMajorityRuleParent(p);
//                }
//                if (mrParent == null || p.getSize() < mrParent.getSize()) {
//                    mrParent = p;
//                }
//            }
//        }
//
//        assert clade.getCredibility() >= 0.5;
//        assert mrParent.getCredibility() >= 0.5;
//
//        clade.setMajorityRuleParent(mrParent);
//        mrParent.addChild(clade);
//
//        return mrParent;
//    }

    private static BiClade findMajorityRuleParent(BiClade clade) {
        if (clade.getMajorityRuleParent() != null) {
            return clade.getMajorityRuleParent();
        }

        Set<BiClade> mrParents = new HashSet<>();

        // find the smallest parent clade that is majority rule
        for (BiClade parent : clade.getParentClades()) {
            if (parent.getCredibility() < 0.5) {
                mrParents.add(findMajorityRuleParent(parent));
            } else {
                mrParents.add(parent);
            }
        }

        BiClade mrParent = null;
        for (BiClade parent : mrParents) {
            if (mrParent == null || parent.getSize() < mrParent.getSize()) {
                mrParent = parent;
            }
        }

        clade.setMajorityRuleParent(mrParent);

        return mrParent;
    }
//    private Set<BiClade> findMajorityRuleConsensusTree(BiClade clade) {
//
//        Set<BiClade> subClades = new HashSet<>();
//
//        for (Pair<BiClade, BiClade> subClade : clade.getSubClades()) {
//            BiClade left = subClade.first;
//            BiClade right = subClade.second;
//            if (left.getCredibility() >= 0.5) {
//                subClades.add(left);
//            } else {
//                subClades.addAll(findMajorityRuleConsensusTree(left));
//            }
//            if (right.getCredibility() >= 0.5) {
//                subClades.add(right);
//            } else {
//                subClades.addAll(findMajorityRuleConsensusTree(right));
//            }
//        }
//
//        if (clade.getCredibility() >= 0.5) {
//            for (BiClade subClade : subClades) {
//                clade.addChild(subClade);
//            }
//            return Set.of(clade);
//        }
//
//        return subClades;
//    }

    private FlexibleNode buildMajorityRuleConsensusTree(BiClade clade) {
        FlexibleNode newNode = new FlexibleNode();
        if (clade.getSize() == 1) {
            newNode.setTaxon(clade.getTaxon());
            newNode.setNumber(clade.getIndex());
        } else {
            for (BiClade subclade: clade.getChildClades()) {
                newNode.addChild(buildMajorityRuleConsensusTree(subclade));
                newNode.setNumber(-1);
            }
        }
        return newNode;
    }
}
