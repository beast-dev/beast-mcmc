/*
 * CharacterTreeKeepFlip.java
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

package dr.evomodel.operators;

import dr.evomodel.substmodel.CharacterTreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * Collapses or restores one non-root node of a sampled character tree, the reversible jump between trees with q and
 * q + 1 clades carried out in fixed dimension. A collapsed node keeps an auxiliary raw rate with the slab density f
 * (the prior on the non-root rates). Restoring draws the rate from h = w Exp(mean s m) + (1 - w) f, with m the slab
 * mean, w = nearZeroWeight and s = nearZeroScale, so half of the proposals (by default) start the node almost
 * transparent; collapsing draws a new auxiliary value from f. The returned log Hastings ratio is the proposal part,
 * log f(a) - log h(b) for a restore and log h(b) - log f(a') for a collapse; the prior on the rates supplies f(b)/f(a).
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeKeepFlip extends SimpleMCMCOperator {

    public static final String CHARACTER_TREE_KEEP_FLIP = "characterTreeKeepFlip";

    private final CharacterTreeModel characterTree;
    private final Parameter rates;
    private final Parameter keep;
    private final ParametricDistributionModel slab;
    private final double nearZeroWeight;
    private final double nearZeroScale;
    private final int[] slots;

    public CharacterTreeKeepFlip(CharacterTreeModel characterTree, Parameter rates, ParametricDistributionModel slab,
                                 double nearZeroWeight, double nearZeroScale, double weight) {
        if (characterTree.getKeep() == null) {
            throw new IllegalArgumentException("the character tree " + characterTree.getId() + " has no keep parameter");
        }
        final int internal = characterTree.getSkeleton().getInternalNodeCount();
        if (rates.getDimension() != internal) {
            throw new IllegalArgumentException("rates has dimension " + rates.getDimension() + "; the character tree has "
                    + internal + " internal nodes");
        }
        if (!(nearZeroWeight >= 0.0 && nearZeroWeight <= 1.0) || !(nearZeroScale > 0.0)) {
            throw new IllegalArgumentException("nearZeroWeight must be in [0, 1] and nearZeroScale positive");
        }
        this.characterTree = characterTree;
        this.rates = rates;
        this.keep = characterTree.getKeep();
        this.slab = slab;
        this.nearZeroWeight = nearZeroWeight;
        this.nearZeroScale = nearZeroScale;
        final int tips = characterTree.getSkeleton().getExternalNodeCount();
        final int rootSlot = characterTree.getSkeleton().getRoot().getNumber() - tips;
        this.slots = new int[internal - 1];
        int k = 0;
        for (int i = 0; i < internal; i++) {
            if (i != rootSlot) {
                slots[k++] = i;
            }
        }
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return CHARACTER_TREE_KEEP_FLIP + "(" + characterTree.getId() + ")";
    }

    private double logProposal(double x) {
        final double nearZeroMean = nearZeroScale * slab.mean();
        return Math.log(nearZeroWeight * Math.exp(-x / nearZeroMean) / nearZeroMean + (1.0 - nearZeroWeight) * slab.pdf(x));
    }

    @Override
    public double doOperation() {
        final int i = slots[MathUtils.nextInt(slots.length)];
        final double old = rates.getParameterValue(i);
        if (keep.getParameterValue(i) == 0.0) {
            final double b = MathUtils.nextDouble() < nearZeroWeight
                    ? MathUtils.nextExponential(1.0 / (nearZeroScale * slab.mean()))
                    : slab.quantile(MathUtils.nextDouble());
            rates.setParameterValue(i, b);
            keep.setParameterValue(i, 1.0);
            return slab.logPdf(old) - logProposal(b);
        }
        final double a = slab.quantile(MathUtils.nextDouble());
        keep.setParameterValue(i, 0.0);
        rates.setParameterValue(i, a);
        return logProposal(old) - slab.logPdf(a);
    }
}
