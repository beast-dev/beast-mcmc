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

import dr.evomodel.branchratemodel.LatentStateBranchRateModel;
import dr.inference.markovjumps.MarkovReward;
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */

public class TwoStateOccupancyMarkovRewardsTest extends MathTestCase {

    private static final double tolerance = 10E-3;

    private double sum(double[] v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x;
        }
        return sum;
    }

    public void testNew() {

        // Equal rates
        double rate = 2.1;
        double prop = 0.5;
        double eps = 0.3;
        double branchLength = 2.4;

        MarkovReward markovReward1 = createMarkovReward(rate, prop);
        MarkovReward markovReward2 = createSericolaMarkovReward(rate, prop);

        double r1 = markovReward1.computePdf(eps * branchLength, branchLength, 0, 0);
        double r2 = markovReward2.computePdf(eps * branchLength, branchLength, 0, 0);

        assertEquals(r1, r2, tolerance);

        // Unequal rates
        prop = 0.501;
        MarkovReward markovReward3 = createMarkovReward(rate, prop);
        MarkovReward markovReward4 = createSericolaMarkovReward(rate, prop);

        double r3 = markovReward3.computePdf(eps * branchLength, branchLength, 0, 0);
        double r4 = markovReward4.computePdf(eps * branchLength, branchLength, 0, 0);

        assertEquals(r3, r4, tolerance);



        System.exit(-1);
    }

    public void testTwoStateSericolaRewards1() {
//        final double rate = 0.0015;
//        final double prop = 0.5;
//        final double eps = 0.01;
//        final double branchLength = 2000.0;

        final double rate = 1;
        final double prop = 0.5;
        final double eps = 0.1;
        final double branchLength = 1.2;

        final boolean print = false;

//
//        TwoStateOccupancyMarkovReward two = (TwoStateOccupancyMarkovReward) markovReward;
//
////        run(markovReward, rate, prop, branchLength, print, 1000);
//
//        System.err.println(markovReward.computePdf(0.5 * branchLength, branchLength, 0, 0));
//
//        System.err.println(new Vector(two.getJumpProbabilities()));
//
//        MarkovReward markovReward2 = createMarkovReward(rate, prop + eps);
//        TwoStateOccupancyMarkovReward two2 = (TwoStateOccupancyMarkovReward) markovReward2;
//
//        System.err.println(markovReward2.computePdf(0.5 * branchLength, branchLength, 0, 0));
//        System.err.println(new Vector(two2.getJumpProbabilities()));
//        System.err.println("");

//        double[][] C = two2.getC();
//        double[][] D = two2.getD();
//
//        System.err.println("C:\n" + new Matrix(C));
//
//        System.err.println("D:\n" + new Matrix(D));

    }

//    public void testTwoStateSericolaRewards1() {
//        final double rate = 0.0015;
////        final double prop = 0.5;
//        final double prop = 0.66666;
//
//        final double branchLength = 2000.0;
//        final boolean print = false;
//
////        MarkovReward markovReward = createMarkovReward(rate, prop);
//        MarkovReward markovReward = createSericolaMarkovReward(rate, prop);
//
//        run(markovReward, rate, prop, branchLength, print, 1000);
//    }

//    public void testTwoStateSericolaRewards2() {
//        final double rate = 0.0015;
//        final double prop = 0.5;
////        final double prop = 0.66666;
//        final double branchLength = 1000.0;
//        final boolean print = false;
//
//        MarkovReward markovReward = createMarkovReward(rate, prop);
////        MarkovReward markovReward = createSericolaMarkovReward(rate, prop);
//
//        run(markovReward, rate, prop, branchLength, print, 1000);
//    }

//    public void testLatentStateBranchRateModel() throws FunctionEvaluationException, MaxIterationsExceededException {
//
//        LatentStateBranchRateModel model = new LatentStateBranchRateModel(
//                new Parameter.Default(0.001), new Parameter.Default(0.5));
//
//        TrapezoidIntegrator integator = new TrapezoidIntegrator();
//
//        final double branchLength = 2000;
//        double integral = integator.integrate(new LatentStateDensityFunction(model, branchLength), 0.0, 1.0);
//
//        System.out.println("testLatentStateBeanchRateModel");
//        System.out.println("Integral = " + integral);
//
//        assertEquals(integral, 1.0, tolerance);
//    }

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

        // TODO Normalization is missing in LatentBranchRateModel
        System.out.println("branchLength = " + branchLength);
        System.out.println("rate = " + rate);
        System.out.println("prop = " + prop);
        System.out.println("Integral = " + sum);
        System.out.println("Mode = " + String.format("%3.2e", modeY) + " at " + modeX);

        assertEquals(sum, 1.0, tolerance);

        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        double integral = 0.0;
        try {
            integral = integrator.integrate(new UnitDensityFunction(markovReward, branchLength, rate, prop), 0.0, 1.0);
        } catch (MaxIterationsExceededException e) {
            e.printStackTrace();
        } catch (FunctionEvaluationException e) {
            e.printStackTrace();
        }
        System.out.println("unt int = " + integral);
        assertEquals(integral, 1.0, tolerance);

        System.out.println("\n");
    }

    private class LatentStateDensityFunction implements UnivariateRealFunction {

        private final LatentStateBranchRateModel model;
        private final double branchLength;

        LatentStateDensityFunction(LatentStateBranchRateModel model, double branchLength) {
            this.model = model;
            this.branchLength = branchLength;
        }

        public double value(double prop) {
            return model.getBranchRewardDensity(prop, branchLength);
        }
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
            return markovReward.computePdf(v, branchLength, 0, 0) /
                    (markovReward.computeConditionalProbability(branchLength, 0, 0)
                            - Math.exp(-rate * prop * branchLength));
        }
    }

    private class UnitDensityFunction implements UnivariateRealFunction {

        private final MarkovReward markovReward;
        private final double branchLength;
        private final double rate;
        private final double prop;

        UnitDensityFunction(MarkovReward markovReward, double branchLength, double rate, double prop) {
            this.markovReward = markovReward;
            this.branchLength = branchLength;
            this.rate = rate;
            this.prop = prop;
        }

        @Override
        public double value(double v) { //throws FunctionEvaluationException {
            double density = markovReward.computePdf(v * branchLength, branchLength, 0, 0) /
                    (markovReward.computeConditionalProbability(branchLength, 0, 0)
                            - Math.exp(-rate * prop * branchLength));
            return density * branchLength;
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
