/*
 * TwoStateOccupancyMarkovRewardsTest.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.substmodel;

import dr.inference.markovjumps.MarkovReward;
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */

public class TwoStateOccupancyMarkovRewardsTest extends MathTestCase {

    public void testTwoStateSericolaRewards1() {
        final double rate = 0.0015;
//        final double prop = 0.5;
        final double prop = 0.66666;

        final double branchLength = 2000.0;
        final boolean print = false;

//        MarkovReward markovReward = createMarkovReward(rate, prop);
        MarkovReward markovReward = createSericolaMarkovReward(rate, prop);

        run(markovReward, rate, prop, branchLength, print, 1000);
    }

    public void testTwoStateSericolaRewards2() {
        final double rate = 0.0015;
//        final double prop = 0.5;
        final double prop = 0.66666;
        final double branchLength = 1000.0;
        final boolean print = false;

//        MarkovReward markovReward = createMarkovReward(rate, prop);
        MarkovReward markovReward = createSericolaMarkovReward(rate, prop);

        run(markovReward, rate, prop, branchLength, print, 1000);
    }

    private void run(MarkovReward markovReward, double rate, double prop, double branchLength,
                     boolean print, int length) {
        DensityFunction densityFunction = new DensityFunction(markovReward, branchLength, rate, prop);

        final double step = branchLength / length;
        int i = 0;
        double sum = 0.0;
        double modeY = 0.0;
        double modeX = 0.0;

        for (double x = 0.0; x <= branchLength; x += step, ++i) {

            double density = 0;
            density = densityFunction.value(x);

            if (x == 0.0) {
                modeY = density;
            } else {
                if (density > modeY) {
                    modeY = density;
                    modeX = x;
                }
            }

            if (x == 0.0 || x == branchLength) {
                sum += density;
            } else {
                sum += 2.0 * density;
            }

            if (print) {
                System.out.println(i + "\t" + String.format("%3.2f", x) + "\t" + String.format("%5.3e", density));
            }
        }
        sum *= (branchLength / 2.0 / length);

        double norm = 1.0 - markovReward.computeCdf(0, branchLength, 0, 0) /
                markovReward.computeConditionalProbability(branchLength, 0, 0);
        sum /= norm; // TODO Normalization is missing in LatentBranchRateModel

        System.out.println("branchLength = " + branchLength);
        System.out.println("rate = " + rate);
        System.out.println("prop = " + prop);
        System.out.println("Integral = " + sum);
        System.out.println("Mode = " + String.format("%3.2e", modeY) + " at " + modeX);
        System.out.println("\n");
    }

    private class DensityFunction implements UnivariateRealFunction {

        private final MarkovReward markovReward;
        private final double branchLength;
        private final double rate;
        private final double prop;

        DensityFunction(MarkovReward markovReward, double branchLength, double rate, double prop) {
            this.markovReward = markovReward;
            this.branchLength = branchLength;
            this.rate = rate;
            this.prop = prop;
        }

        @Override
        public double value(double v) { //throws FunctionEvaluationException {
            return markovReward.computePdf(v, branchLength, 0, 0)
                    / markovReward.computeConditionalProbability(branchLength, 0, 0);
        }
    }

    private double[] createLatentInfinitesimalMatrix(final double rate, final double prop) {
        double[] mat = new double[]{
                -rate * prop, rate * prop,
                rate * (1.0 - prop), -rate * (1.0 - prop)
        };
        return mat;
    }

    private SericolaSeriesMarkovReward createSericolaMarkovReward(final double rate, final double prop) {
        double[] r = new double[]{0.0, 1.0};
        return new SericolaSeriesMarkovReward(createLatentInfinitesimalMatrix(rate, prop), r, 2);
    }

    private TwoStateOccupancyMarkovReward createMarkovReward(final double rate, final double prop) {
        TwoStateOccupancyMarkovReward markovReward = new
                TwoStateOccupancyMarkovReward(
                createLatentInfinitesimalMatrix(rate, prop)
        );
        return markovReward;
    }
}
