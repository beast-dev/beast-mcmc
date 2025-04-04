/*
 * TreeSummaryStatistic.java
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

package dr.app.treestat.statistics;

import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

import java.util.Map;

/**
 * An interface and collection of tree summary statistics.
 *
 *
 * @author Alexei Drummond
 */
public interface TreeSummaryStatistic extends SummaryStatisticDescription {

	int getStatisticDimensions(Tree tree);

	String getStatisticLabel(Tree tree, int i);

    void setTaxonList(TaxonList taxonList);
    void setInteger(int value);
    void setDouble(double value);
    void setString(String value);

    /**
	 * @return the value of this summary statistic for the given tree.
	 */
	double[] getSummaryStatistic(Tree tree);

	public abstract class Factory implements SummaryStatisticDescription {
		public TreeSummaryStatistic createStatistic() {
			throw new RuntimeException("This factory method is not implemented");
		}

		public boolean allowsWholeTree() { return true; }
		public boolean allowsTaxonList() { return false; }
		public boolean allowsInteger() { return false; }
		public boolean allowsDouble() { return false; }
        public boolean allowsString() { return false; }

		public String getValueName() { return ""; }

    }
}

