/*
 * BirthDeathCollapseNClustersStatistic.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/*
        This file is part of BEAST.

        BEAST is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        BEAST is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with DISSECT.  If not, see <http://www.gnu.org/licenses/>.
*/

package dr.evomodel.alloppnet.speciation;


import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodel.alloppnet.parsers.BirthDeathCollapseNClustersStatisticParser;
import dr.inference.model.Statistic;

/**
 * @author Graham Jones
 *         Date: 01/10/2013
 */
public class BirthDeathCollapseNClustersStatistic extends Statistic.Abstract {
    private SpeciesTreeModel spptree;
    private BirthDeathCollapseModel bdcm;


    public BirthDeathCollapseNClustersStatistic(SpeciesTreeModel spptree, BirthDeathCollapseModel bdcm) {
        super(BirthDeathCollapseNClustersStatisticParser.BDC_NCLUSTERS_STATISTIC);
        this.spptree = spptree;
        this.bdcm = bdcm;
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        int ninodes = spptree.getInternalNodeCount();
        int n =  0;
        for (int i = 0; i < ninodes; i++) {
            double h = spptree.getNodeHeight(spptree.getInternalNode(i));
            if (!BirthDeathCollapseModel.belowCollapseHeight(h, bdcm.getCollapseHeight())) {
                n++;
            }
        }
        return n+1;
    }
}
