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
import dr.util.Pair;

import java.util.*;

public class MajorityRuleTreeBuilder {

    public MutableTree getMajorityRuleConsensusTree(CladeSystem cladeSystem, TaxonList taxonList) {
        BiClade rootClade = (BiClade)cladeSystem.getRootClade();

        Set<Clade> majorityRuleClades = cladeSystem.getTopClades(0.5);
        findMajorityRuleConsensusTree(rootClade, majorityRuleClades);

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

    private BitSet findMajorityRuleConsensusTree(BiClade clade) {

        Set<BiClade> subClades = new HashSet<>();

        for (Pair<BiClade, BiClade> subClade : clade.getSubClades()) {
            if (subClade.first.getSize() > 1) {
                Set<BiClade> leftSubClades = findMajorityRuleConsensusTree(subClade.first, majorityRuleClades);
                if (subClade.first.getCredibility() >= 0.5) {
                    subClades.add(subClade.first);
                } else {
                    // may be no subclades (whi
                    subClades.addAll(leftSubClades);
                }
            }
            Set<BiClade> rightSubClades = findMajorityRuleConsensusTree(subClade.second, majorityRuleClades);
            if (subClade.second.size > 1 && subClade.second.getCredibility() >= 0.5) {
                subClades.add(subClade.second);
            } else {
                subClades.addAll(rightSubClades);
            }
        }

        clade.majorityRuleSubClades = subClades;

        return subClades;
    }

    private FlexibleNode buildMajorityRuleConsensusTree(BiClade clade) {
        FlexibleNode newNode = new FlexibleNode();
        if (clade.getSize() == 1) {
            newNode.setTaxon(clade.getTaxon());
            newNode.setNumber(clade.getIndex());
        } else {
            for (BiClade subclade: clade.majorityRuleSubClades) {
                newNode.addChild(buildMajorityRuleConsensusTree(subclade));
                newNode.setNumber(-1);
            }
        }
        return newNode;
    }
}
