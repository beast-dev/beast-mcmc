/*
 * CoveringSubtreeSizeStatistic.java
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
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * @author Alexei Drummond
 *
 * @version $Id: CoveringSubtreeSizeStatistic.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 */
class CoveringSubtreeSizeStatistic extends AbstractTreeSummaryStatistic {

	private CoveringSubtreeSizeStatistic(Map characterMap, String characterState) {
		this.characterState = characterState;
		Iterator i = characterMap.entrySet().iterator();
		while (i.hasNext()) {
			Map.Entry entry = (Map.Entry)i.next();
			if (entry.getValue().equals(characterState)) {
				Object o = entry.getKey();
				if (o instanceof String) {
					taxa.addTaxon(new Taxon((String)o));
				} else if (o instanceof Taxon) {
					taxa.addTaxon((Taxon)o);
				} else throw new RuntimeException("Unknown key type!");
			}
		}
	}

	public double[] getSummaryStatistic(Tree tree) {
		try {
			Set leafSet = Tree.Utils.getLeavesForTaxa(tree, taxa);
			NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leafSet);
			return new double[] { Tree.Utils.getLeafCount(tree, node) };

		} catch (Tree.MissingTaxonException e) {
			throw new RuntimeException("Missing taxon!");
		}
	}

	public String getSummaryStatisticName() {
		return "coveringSubtreeSize(" + characterState + ")";
	}

	public String getSummaryStatisticDescription() {
		return "The number of taxa in the smallest subtree that completely covers the character state " + characterState;
	}

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }

	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic(Map characterMap, String characterState) {
			return new CoveringSubtreeSizeStatistic(characterMap, characterState);
		}

		public String getSummaryStatisticName() {
			return "coveringSubtreeSize";
		}

		public String getSummaryStatisticDescription() {
			return "The number of taxa in the smallest subtree that completely covers a character state";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.GENERAL; }

		public boolean allowsWholeTree() { return false; };
		public boolean allowsCharacter() { return false; };
		public boolean allowsCharacterState() { return true; };
		public boolean allowsTaxonList() { return false; };
		public boolean allowsInteger() { return false; };
		public boolean allowsDouble() { return false; };
	};

	private String characterState;
	private Taxa taxa = new Taxa();
}