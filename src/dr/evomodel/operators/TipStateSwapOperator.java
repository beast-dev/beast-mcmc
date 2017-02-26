/*
 * TipStateSwapOperator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.operators;

import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Marc A. Suchard
 */
public class TipStateSwapOperator extends SimpleMCMCOperator {

    public static final String TIP_STATE_OPERATOR = "tipStateSwapOperator";

    public TipStateSwapOperator(AncestralStateBeagleTreeLikelihood treeLikelihood, double weight) {
        this.treeLikelihood = treeLikelihood;
        setWeight(weight);
        int patternCount = treeLikelihood.getPatternCount();
        states1 = new int[patternCount];
        states2 = new int[patternCount];
    }

    private final int[] states1;
    private final int[] states2;
    private int index1;
    private int index2;

    public double doOperation() {

        int tipCount = treeLikelihood.getTreeModel().getExternalNodeCount();

        // Choose two tips to swap
        index1 = MathUtils.nextInt(tipCount);
        index2 = index1;
        while (index2 == index1)
            index2 = MathUtils.nextInt(tipCount);

        swap(index1, index2);

        treeLikelihood.makeDirty();

        return 0;
    }

    private void swap(int i, int j) {

        treeLikelihood.getTipStates(i, states1);
        treeLikelihood.getTipStates(j, states2);

        treeLikelihood.setTipStates(j, states1);
        treeLikelihood.setTipStates(i, states2);
    }

    public void reject() {
        super.reject();
        // There is currently no restore functions for tip states, so manually adjust state
        swap(index1, index2);
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String getOperatorName() {
        return TIP_STATE_OPERATOR;
    }

    private final AncestralStateBeagleTreeLikelihood treeLikelihood;
}