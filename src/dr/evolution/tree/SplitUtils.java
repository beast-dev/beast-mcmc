/*
 * SplitUtils.java
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

package dr.evolution.tree;

import dr.evolution.util.TaxonList;

/**
 * utilities for split systems
 *
 * @version $Id: SplitUtils.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Korbinian Strimmer
 */
public class SplitUtils
{
	//
	// Public stuff
	//

	/**
	 * creates a split system from a tree
	 * (using a pre-specified order of sequences)
	 *
	 * @param taxonList the list of taxa (order is important)
	 * @param tree
	 */
	public static SplitSystem getSplits(TaxonList taxonList, Tree tree)
	{
		int size = tree.getInternalNodeCount()-1;
		SplitSystem splitSystem = new SplitSystem(taxonList, size);

		boolean[][] splits = splitSystem.getSplitVector();
	
		int j = 0;
		for (int i = 0; i < tree.getInternalNodeCount(); i++) {
			NodeRef node = tree.getInternalNode(i);
			if (node != tree.getRoot()) {
				getSplit(taxonList, tree, node, splits[j]);
				j++;
			}
		}

		return splitSystem;
	}



	/**
	 * creates a split system from a tree
	 * (using tree-induced order of sequences)
	 *
	 * @param tree
	 */
	public static SplitSystem getSplits(Tree tree)
	{
        return getSplits(tree, tree);
	}



	/**
	 * get split for branch associated with internal node
	 *
	 * @param taxonList
     * @param tree
	 * @param internalNode Node
	 * @param split
	 */
	public static void getSplit(TaxonList taxonList, Tree tree, NodeRef internalNode, boolean[] split)
	{
		if (tree.isExternal(internalNode) || tree.isRoot(internalNode))
		{
			throw new IllegalArgumentException("Only internal nodes (and no root) nodes allowed");
		}

		// make sure split is reset
		for (int i = 0; i < split.length; i++)
		{
			split[i] = false;
		}

		// mark all leafs downstream of the node

		for (int i = 0; i < tree.getChildCount(internalNode); i++)
		{
			markNode(taxonList, tree, internalNode, split);
		}

		// standardize split (i.e. first index is alway true)
		if ( !split[0] )
		{
			for (int i = 0; i < split.length; i++)
			{
                split[i] = !split[i];
			}
		}
	}

	/**
	 * checks whether two splits are identical
	 * (assuming they are of the same length
	 * and use the same leaf order)
	 *
	 * @param s1 split 1
	 * @param s2 split 2
	 */
	public static boolean isSame(boolean[] s1, boolean[] s2)
	{
		boolean reverse;
        reverse = s1[0] != s2[0];

		if (s1.length != s2.length)
			throw new IllegalArgumentException("Splits must be of the same length!");

		for (int i = 0; i < s1.length; i++)
		{
			if (reverse)
			{
				// splits not identical
				if (s1[i] == s2[i]) return false;
			}
			else
			{
				// splits not identical
				if (s1[i] != s2[i]) return false;
			}
		}

		return true;
	}

	//
	// Package stuff
	//

	static void markNode(TaxonList taxonList, Tree tree, NodeRef node, boolean[] split)
	{
		if (tree.isExternal(node))
		{
			String name = tree.getTaxonId(node.getNumber());
			int index = taxonList.getTaxonIndex(name);

			if (index < 0)
			{
				throw new IllegalArgumentException("INCOMPATIBLE IDENTIFIER (" + name + ")");
			}

			split[index] = true;
		}
		else
		{
			for (int i = 0; i < tree.getChildCount(node); i++)
			{
				markNode(taxonList, tree, tree.getChild(node, i), split);
			}
		}
	}

}
