/*
 * ImportanceDistribution.java
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

/**
 * 
 */
package dr.evolution.tree;

/**
 * @author Sebastian Hoehna
 *
 */
public interface ImportanceDistribution {
	
	/**
	 * 
	 * @param tree
	 *            - the tree to be added
	 */
	public void addTree(Tree tree);	
	
	/**
	 * 
	 * Splits a clade into two sub-clades according to the importance distribution
	 * 
	 * @param parent - the clade which is split
	 * @param children - a call by reference parameter which is an empty, two element array of clades at time of call and contains the to sub clades afterwards
	 * @return the chance for this split
	 */
	public double splitClade(Clade parent, Clade[] children);
	
	/**
	 * 
	 * Calculates the probability of a given tree.
	 * 
	 * @param tree
	 *            - the tree to be analyzed
	 * @return estimated posterior probability in log
	 */
	public double getTreeProbability(Tree tree);

}
