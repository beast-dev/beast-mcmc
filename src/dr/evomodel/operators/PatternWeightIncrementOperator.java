/*
 * PatternWeightIncrementOperator.java
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

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import jebl.math.Random;

/**
 * @author Andrew Rambaut
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class PatternWeightIncrementOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String PATTERN_WEIGHT_INCREMENT_OPERATOR = "patternWeightIncrementOperator";

    public PatternWeightIncrementOperator(BeagleTreeLikelihood treeLikelihood, double weight) {
        this.treeLikelihood = treeLikelihood;
        setWeight(weight);

        finalPatternWeights = treeLikelihood.getPatternWeights();

        double[] weights = new double[finalPatternWeights.length];
        treeLikelihood.setPatternWeights(finalPatternWeights);
//        treeLikelihood.setPatternWeights(weights);
    }

    public double doOperation() {

        if (allPatternsAdded) {
            return 0.0;
        }

        double[] weights = treeLikelihood.getPatternWeights();
        double[] w = new double[weights.length];

        double sum = 0.0;
        for (int i = 0; i < weights.length; i++) {
            w[i] = finalPatternWeights[i] - weights[i];
            sum += w[i];
        }

        // System.out.println("PatternWeightIncrementOperator - Sites remaining: " + sum);
        if (sum < 1.0) {
            allPatternsAdded = true;
            System.out.println("PatternWeightIncrementOperator - All sites added");
            return 0.0;
        }

        for (int i = 0; i < weights.length; i++) {
            w[i] /= sum;
            if (i > 0) {
                w[i] += w[i - 1];
            }
        }

        int r = Random.randomChoice(w);

        weights[r] ++;

        if (weights[r] > finalPatternWeights[r]) {
            throw new RuntimeException("Pattern weight exceeding final weight");
        }

        treeLikelihood.setPatternWeights(weights);

        treeLikelihood.makeDirty();

        return 0;
    }

    public void reject() {
        super.reject();
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
        return PATTERN_WEIGHT_INCREMENT_OPERATOR;
    }

    private final BeagleTreeLikelihood treeLikelihood;
    private final double[] finalPatternWeights;

    private boolean allPatternsAdded = false;

}