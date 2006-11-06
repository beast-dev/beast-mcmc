/*
 * CladeSet.java
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

package dr.evolution.tree;

import dr.evolution.util.TaxonList;
import dr.util.FrequencySet;

import java.util.BitSet;

/**
 * Stores a set of unique clades for a tree 
 *
 * @version $Id: CladeSet.java,v 1.8 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class CladeSet extends FrequencySet
{
	//
	// Public stuff
	//

	public CladeSet()
	{
	}

	/**
	 * @param tree
	 */
	public CladeSet(Tree tree)
	{
		this.taxonList = tree;
		add(tree);
	}

	/**
	 * @param taxonList  a set of taxa used to label the tips
	 */
	public CladeSet(Tree tree, TaxonList taxonList)
	{
		this.taxonList = taxonList;
		add(tree);
	}

	/** get number of unique clades */
	public int getCladeCount()
	{		
		return size();
	}

	/** get clade bit set */
	public String getClade(int index)
	{
		BitSet bits = (BitSet)get(index);
		
		StringBuffer buffer = new StringBuffer("{");
		boolean first = true;
		for (int i = 0; i < bits.length(); i++) {
			if (bits.get(i)) {
				if (!first) {
					buffer.append(", ");					
				} else {
					first = false;
				}
				buffer.append(taxonList.getTaxonId(i));
			}
		}
		buffer.append("}");
		return buffer.toString();
	}

	/** get clade frequency */
	public int getCladeFrequency(int index)
	{		
		return getFrequency(index);
	}

	/** adds all the clades in the tree */
	public void add(Tree tree)
	{		
		if (taxonList == null) {
			taxonList = tree;
		}
		
		// Recurse over the tree and add all the clades (or increment their
		// frequency if already present). The root clade is not added.
		addClades(tree, tree.getRoot(), null);
	}

	private void addClades(Tree tree, NodeRef node, BitSet bits) {
	
		if (tree.isExternal(node)) {
		
			if (taxonList != null) {
				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits.set(index);
			} else {
				bits.set(node.getNumber());
			}
		} else {
			
			BitSet bits2 = new BitSet();
			for (int i = 0; i < tree.getChildCount(node); i++) {
			
				NodeRef node1 = tree.getChild(node, i);

				addClades(tree, node1, bits2);
			}
			
			add(bits2, 1);
			
			if (bits != null) {
				bits.or(bits2);
			}
		}
	}
		
	/** adds all the clades in the CladeSet */
	public void add(CladeSet cladeSet)
	{		
		for (int i = 0, n = cladeSet.getCladeCount(); i < n; i++) {
			add(cladeSet.getClade(i), cladeSet.getCladeFrequency(i));
		}
	}
	
	public boolean hasClade(int index, Tree tree) {
		BitSet bits = (BitSet)get(index);

		NodeRef[] mrca = new NodeRef[1];
		findClade(bits, tree, tree.getRoot(), mrca);			
	
		return (mrca[0] != null);
	}
	
	public int findClade(BitSet bitSet, Tree tree, NodeRef node, NodeRef[] cladeMRCA) {
		
		if (tree.isExternal(node)) {
		
			if (taxonList != null) {
				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				if (bitSet.get(index)) return 1;
			} else {
				if (bitSet.get(node.getNumber())) return 1;
			} 
			return -1;
		} else {
			int count = 0;
			for (int i = 0; i < tree.getChildCount(node); i++) {
			
				NodeRef node1 = tree.getChild(node, i);

				int childCount = findClade(bitSet, tree, node1, cladeMRCA);
				
				if (childCount != -1 && count != -1) {
					count += childCount;
				} else count = -1;
			}
			
			if (count == bitSet.cardinality()) cladeMRCA[0] = node;
			
			return count;
		}
	}

	//
	// Private stuff
	//
	TaxonList taxonList = null;
}
