/*
 * TMRCAStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.Statistic;

import java.util.Set;

/**
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 */
public class AgeStatistic extends Statistic.Abstract {

    public AgeStatistic(String name, Statistic heightStatistic) {
        super(name);
        this.heightStatistic = heightStatistic;
        if (Taxon.getMostRecentDate() != null) {
            isBackwards = Taxon.getMostRecentDate().isBackwards();
            mostRecentTipTime = Taxon.getMostRecentDate().getAbsoluteTimeValue();
        } else {
            // give node heights or taxa don't have dates
            mostRecentTipTime = Double.NaN;
            isBackwards = false;
        }
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {

        if (!Double.isNaN(mostRecentTipTime)) {
            if (isBackwards) {
                return mostRecentTipTime + heightStatistic.getStatisticValue(dim);
            } else {
                return mostRecentTipTime - heightStatistic.getStatisticValue(dim);
            }
        } else {
            return Double.NaN;
        }
    }

    private Statistic heightStatistic = null;
    private final double mostRecentTipTime;
    private final boolean isBackwards;

}
