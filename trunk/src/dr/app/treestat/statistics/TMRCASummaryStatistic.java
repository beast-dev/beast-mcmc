/*
 * TMRCASummaryStatistic.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.treestat.statistics;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.*;

import java.util.*;

/**
 *
 * @version $Id: TMRCASummaryStatistic.java,v 1.3 2006/05/09 10:24:27 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class TMRCASummaryStatistic extends AbstractTreeSummaryStatistic {

	private TMRCASummaryStatistic(Map<String, Object> characterMap, String characterState) {
		taxonList = new Taxa();
        this.characterState = characterState;
        for (Map.Entry entry : characterMap.entrySet()) {
            if (entry.getValue().equals(characterState)) {
                Object o = entry.getKey();
                if (o instanceof String) {
                    ((Taxa) taxonList).addTaxon(new Taxon((String) o));
                } else if (o instanceof Taxon) {
                    ((Taxa) taxonList).addTaxon((Taxon) o);
                } else throw new RuntimeException("Unknown key type!");
            }
        }
    }

	private TMRCASummaryStatistic(TaxonList taxonList) {
		this.taxonList = taxonList;
	}

	private TMRCASummaryStatistic() {
	}

	public double[] getSummaryStatistic(Tree tree) {
		if (taxonList == null) {
            return new double[] { tree.getNodeHeight(tree.getRoot()) };
        }
		try {
			Set<String> leafSet = Tree.Utils.getLeavesForTaxa(tree, taxonList);
			NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leafSet);
			if (node == null) throw new RuntimeException("No node found that is MRCA of " + leafSet);
			return new double[] { tree.getNodeHeight(node) };
		} catch (Tree.MissingTaxonException e) {
			throw new RuntimeException("Missing taxon!");
		}
	}

	public String getSummaryStatisticName() {
		if (characterState != null) {
			return "tMRCA(" + characterState + ")";
		} else if (taxonList != null) {
			return "tMRCA(" + taxonList.getId() + ")";
		} else {
			return "tMRCA";
		}
	}

	public String getSummaryStatisticDescription() {
		if (characterState != null) {
			return "The time of the most recent common ancestor of the character state " + characterState;
		} else if (taxonList != null) {
			return "The time of the most recent common ancestor of the given taxon list";
		}
		return "The time of the most recent common ancestor of a set of taxa. In order to use this statistic, a taxon set must be defined (see the Taxon Set tab).";
	}

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }

	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic(Map<String, Object> characterMap, String characterState) {
			return new TMRCASummaryStatistic(characterMap, characterState);
		}

		public TreeSummaryStatistic createStatistic(TaxonList taxonList) {
			return new TMRCASummaryStatistic(taxonList);
		}

		public TreeSummaryStatistic createStatistic() {
			return new TMRCASummaryStatistic();
		}

		public String getSummaryStatisticName() {
			return "tMRCA";
		}

		public String getSummaryStatisticDescription() {
			return "The time of the most recent common ancestor";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.GENERAL; }

		public boolean allowsWholeTree() { return true; }

        public boolean allowsCharacter() { return false; }

        public boolean allowsCharacterState() { return false; }

        public boolean allowsTaxonList() { return true; }
    };

	private String characterState = null;
	private TaxonList taxonList = null;
}

