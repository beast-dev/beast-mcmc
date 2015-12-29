/*
 * SummaryStatisticDescription.java
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

package dr.app.treestat.statistics;

/**
 * An interface and collection of tree summary statistics.
 *
 * @version $Id: SummaryStatisticDescription.java,v 1.1 2005/09/26 22:14:15 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public interface SummaryStatisticDescription {

	String getSummaryStatisticName();

	String getSummaryStatisticDescription();

	String getSummaryStatisticReference();

	boolean allowsPolytomies();

	boolean allowsNonultrametricTrees();

	boolean allowsUnrootedTrees();

	Category getCategory();

	public class Category {

		public static final Category TREE_SHAPE = new Category("Tree shape");
		public static final Category PHYLOGENETIC = new Category("Phylogenetic");
		public static final Category POPULATION_GENETIC = new Category("Population genetic");
		public static final Category GENERAL = new Category("General");
		public static final Category SPECIATION = new Category("Speciation/Birth-death");

		private Category(String name) {
			this.name = name;
		}

		public String toString() { return name; }

		private String name;
	}

}

