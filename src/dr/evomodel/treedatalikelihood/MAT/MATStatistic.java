/*
 * MATStatistic.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.inference.model.Statistic;

/**
 * Logs the number of mutations and ambiguities in the current MAT, counted per site (i.e. weighted
 * by pattern weight).
 */
public class MATStatistic extends Statistic.Abstract {

    private static final String[] DIMENSION_NAMES = {"mutations", "ambiguities"};

    private final MATDataLikelihoodDelegate delegate;

    public MATStatistic(String name, MATDataLikelihoodDelegate delegate) {
        super(name);
        this.delegate = delegate;
    }

    @Override
    public int getDimension() {
        return DIMENSION_NAMES.length;
    }

    @Override
    public String getDimensionName(int dim) {
        return (getStatisticName() == null ? "" : getStatisticName() + ".") + DIMENSION_NAMES[dim];
    }

    @Override
    public double getStatisticValue(int dim) {
        MutationStore store = delegate.getMutationStore();
        double[] weights = delegate.getPatternWeights();
        double count = 0.0;
        for (int n = 0; n < store.getNodeCount(); n++) {
            BranchMutations mutations = store.getBranchMutations(n);
            for (int k = 0; k < mutations.size(); k++) {
                if (mutations.isAmbiguity(k) == (dim == 1)) {
                    count += weights[mutations.getPattern(k)];
                }
            }
        }
        return count;
    }
}
