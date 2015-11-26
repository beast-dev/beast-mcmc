/*
 * MultiEpochExponentialTest.java
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

package test.dr.evomodel.coalescent;

import dr.evolution.coalescent.ExponentialExponential;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.coalescent.MultiEpochExponential;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evomodel.coalescent.VariableDemographicModel;
import junit.framework.TestCase;

/**
 * @author Marc A. Suchard
 */
public class MultiEpochExponentialTest extends TestCase {

    public void testExponentialExponential() {

        Units.Type units = Units.Type.YEARS;

        ExponentialExponential ee = new ExponentialExponential(units);
        ee.setN0(N0);
        ee.setGrowthRate(rates[0]);
        ee.setAncestralGrowthRate(rates[1]);
        ee.setTransitionTime(transitionTimes[0]);

        MultiEpochExponential mee = new MultiEpochExponential(units, 2);
        mee.setN0(N0);
        for (int i = 0;i < rates.length; ++i) {
            mee.setGrowthRate(i, rates[i]);
        }
        for (int i = 0; i < transitionTimes.length; ++i) {
            mee.setTransitionTime(i, transitionTimes[i]);
        }

        for (double time = 0; time < 20; time += 1.0) {
            double eeDemo = ee.getDemographic(time);
            double meeDemo = mee.getDemographic(time);
            assertEquals(eeDemo, meeDemo, tolerance1);
        }

        double start = 0.0;
        double finish = 1.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = ee.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }


        start = 0.5;
        finish = 1.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = ee.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }

        start = 11.0;
        finish = 11.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = ee.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }
    }


    public void testThreeExponential() {

        Units.Type units = Units.Type.YEARS;

        ExponentialGrowth e = new ExponentialGrowth(units);
        e.setN0(N0);
        e.setGrowthRate(rates3[0]);

        MultiEpochExponential mee = new MultiEpochExponential(units, 3);
        mee.setN0(N0);
        for (int i = 0; i < rates3.length; ++i) {
            mee.setGrowthRate(i, rates3[i]);
        }
        for (int i = 0; i < transitionTimes3.length; ++i) {
            mee.setTransitionTime(i, transitionTimes3[i]);
        }

        double start = 0.0;
        double finish = 1.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = e.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }


        start = 0.5;
        finish = 1.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = e.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }

        start = 11.0;
        finish = 11.0;

        for (; finish < 20.0; finish += 1.0) {
            double eeInt = e.getIntegral(start, finish);
            double meeIntN = mee.getNumericalIntegral(start, finish);
            double meeIntA = mee.getAnalyticIntegral(start, finish);
//            System.err.println(finish + ": " + eeInt + " " + meeIntN + " " + meeIntA);
            assertEquals(eeInt, meeIntN, tolerance1);
            assertEquals(meeIntN, meeIntA, tolerance2);
        }
    }


    double N0 = 100;

    double[] rates = new double[] { 0.2, -0.2 };
    double[] transitionTimes = new double[] { 10.0 };

    double[] rates3 = new double[] { 0.1, 0.1, 0.1 };
    double[] transitionTimes3 = new double[] { 10.0, 15.0 };

    private final static double tolerance1 = 1E-10;
    private final static double tolerance2 = 1E-5;
}
