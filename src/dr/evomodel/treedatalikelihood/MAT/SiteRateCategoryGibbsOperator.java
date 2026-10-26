/*
 * SiteRateCategoryGibbsOperator.java
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

import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;

/**
 * A Gibbs update of the site rate categories of a MAT likelihood with a discrete site rate model: every
 * variable site's category is drawn from its full conditional given its history and the tree (see
 * {@link MATDataLikelihoodDelegate#resampleCategories}).
 * <p>
 * With more than one rate category the MAT samples each variable site's category, as it samples its
 * history, so that its likelihood still decomposes by branch; summed over the categories it is that of
 * the site rate model. The categories only change through this operator.
 * <p>
 * Each update draws every site's category and makes the whole likelihood be recomputed, and the categories
 * mix quickly, so a small weight is enough: about 1% of the subtree moves' weight (on small-2, 10% of it
 * took nearly half the run time for no gain).
 */
public class SiteRateCategoryGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String OPERATOR_NAME = "siteRateCategoryGibbs";

    private final MATDataLikelihoodDelegate delegate;

    public SiteRateCategoryGibbsOperator(MATDataLikelihoodDelegate delegate, double weight) {
        if (delegate.getCategoryCount() <= 1) {
            throw new IllegalArgumentException("The MAT likelihood has no site rate categories to sample: its site " +
                    "model has a single rate category");
        }
        this.delegate = delegate;
        delegate.setCategoryOperator();
        setWeight(weight);
    }

    @Override
    public double doOperation() {
        delegate.resampleCategories();
        return 0.0;
    }

    @Override
    protected String getBaseOperatorName() {
        return OPERATOR_NAME + "(" + delegate.getTree().getId() + ")";
    }
}
