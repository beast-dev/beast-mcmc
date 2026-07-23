/*
 * RobinsonFouldsMetric.java
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

package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Clade;
import dr.evolution.tree.Tree;

import java.util.*;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Andrew Rambaut
 */
public class RobinsonFouldsMetric implements TreeMetric {
	public static Type TYPE = Type.ROBINSON_FOULDS;

	public RobinsonFouldsMetric() {

	}

	@Override
	public double getMetric(Tree tree1, Tree tree2) {

		checkTreeTaxa(tree1, tree2);

		Set<Clade> clades1 = Clade.getCladeSet(tree1);
		Set<Clade> clades2 = Clade.getCladeSet(tree2);

		clades1.removeAll(clades2);

		// Technically RF would be twice this because it doesn't assume
		// the same set of tips in both trees (so may have a different
		// number of clades missing from each).
		return clades1.size();
	}

	@Override
	public Type getType() {
		return TYPE;
	}

	@Override
	public String toString() {
		return getType().getShortName();
	}

	// todo - add in Citable:
	// Robinson, D. R.; Foulds, L. R. (1981). "Comparison of phylogenetic trees". Mathematical Biosciences. 53: 131-147. doi:10.1016/0025-5564(81)90043-2.
}
