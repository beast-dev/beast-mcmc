/*
 * TipStateOperator.java
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
 */

package dr.evomodel.operators;

import dr.evolution.util.Taxon;
import dr.evomodel.tipstatesmodel.TimeVaryingFrequenciesModel;
import dr.evomodel.treedatalikelihood.TipStateAccessor;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.List;

import static dr.evomodelxml.operators.TipStateOperatorParser.TIP_STATE_OPERATOR;

/**
 * @author Marc A. Suchard
 */
public class TipStateOperator extends SimpleMCMCOperator {

    private final Taxon taxon;
    private final List<TipStateAccessor> treeLikelihoods;
    private final TimeVaryingFrequenciesModel frequencies;
    private final int patternCount;

    private final static boolean ONLy_USE_WITH_TIME_VARYING_MODEL = false;

    private int[] previousStates;

    private static final boolean DEBUG = false;

    public TipStateOperator(TimeVaryingFrequenciesModel frequencies,
                            double weight) {

        this.treeLikelihoods = frequencies.getAccessors();
        this.taxon = frequencies.getTaxon();
        this.frequencies = frequencies;
        this.patternCount = treeLikelihoods.get(0).getPatternCount();

        setWeight(weight);
    }

    private void setCurrentStates(int tipNum, int[] states) {
        for (TipStateAccessor accessor : treeLikelihoods) {
            accessor.setTipStates(tipNum, states);
        }
    }

    private int[] getCurrentStates(int tipNum) {
        int[] currentStates = new int[patternCount];
        treeLikelihoods.get(0).getTipStates(tipNum, currentStates);

        if (treeLikelihoods.size() > 1) {
            int[] test = new int[patternCount];
            for (int i = 1; i < treeLikelihoods.size(); ++i) {
                treeLikelihoods.get(i).getTipStates(tipNum, test);
                if (notEqual(currentStates, test)) {
                   throw new RuntimeException("Inconsistent states");
                }
            }
        }

        return currentStates;
    }

    private static boolean notEqual(int[] lhs, int[] rhs) {
        if (lhs.length != rhs.length) {
            return true;
        }
        for (int i = 0; i < lhs.length; ++i) {
            if (lhs[i] != rhs[i]) {
                return true;
            }
        }
        return false;
    }

    private int[] sample(double[] probabilities) {
        int[] newStates = new int[patternCount];
        for (int i = 0; i < patternCount; ++i) {
            newStates[i] = MathUtils.randomChoicePDF(probabilities);
        }
        return newStates;
    }

    public double doOperation() {

        double[] probabilities = frequencies.getProbabilities(taxon);
        int taxonIndex = frequencies.getTipIndex(taxon);

        previousStates = getCurrentStates(taxonIndex);
        int[] newStates = sample(probabilities);
        setCurrentStates(taxonIndex, newStates);

        if (notEqual(previousStates, newStates)) {
            frequencies.makeDirty();

            if (DEBUG) {
                System.err.println("state change");
            }
        }

        double logRatio = 0;
        for (int i = 0; i < patternCount; ++i) {
            if (DEBUG) {
                System.err.println("\t " + previousStates[i] + " -> " + newStates[i]);
            }

            if (previousStates[i] < probabilities.length) {
                logRatio += Math.log(probabilities[previousStates[i]]);
            }
            logRatio -= Math.log(probabilities[newStates[i]]);
        }

        return ONLy_USE_WITH_TIME_VARYING_MODEL ? 0.0 : logRatio;
    }

    @Override
    public String getOperatorName() {
        return TIP_STATE_OPERATOR;
    }

    public void accept(double deviation) {
        super.accept(deviation);

        if (DEBUG) {
            System.err.println("accept");
        }
    }

    public void reject() {
        super.reject();

        if (DEBUG) {
            System.err.println("reject");
        }

        int taxonIndex = frequencies.getTipIndex(taxon);
        setCurrentStates(taxonIndex, previousStates);
    }
}