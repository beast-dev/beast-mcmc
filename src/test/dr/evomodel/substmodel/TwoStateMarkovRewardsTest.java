/*
 * TwoStateMarkovRewardsTest.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */

public class TwoStateMarkovRewardsTest extends MathTestCase {

    public void testTwoStateRewards() {
        DataType dataType = TwoStates.INSTANCE;
        FrequencyModel freqModel = new FrequencyModel(TwoStates.INSTANCE, new double[]{0.5, 0.5});
        Parameter rates = new Parameter.Default(new double[]{4.0, 6.0});
        ComplexSubstitutionModel twoStateModel = new ComplexSubstitutionModel("two", dataType, freqModel, rates) {
//    protected EigenSystem getDefaultEigenSystem(int stateCount) {
//        return new DefaultEigenSystem(stateCount);
//    }

        };
        twoStateModel.setNormalization(false);

        MarkovJumpsSubstitutionModel markovRewards = new MarkovJumpsSubstitutionModel(twoStateModel,
                MarkovJumpsType.REWARDS);

        double[] r = new double[2];
        double[] q = new double[4];
        double[] c = new double[4];

        int mark = 0;
        double weight = 1.0;

        r[mark] = weight;
        markovRewards.setRegistration(r);

        twoStateModel.getInfinitesimalMatrix(q);
        System.out.println("Q = " + new Vector(q));

        System.out.println("Reward for state 0");
        double time = 1.0;
        markovRewards.computeCondStatMarkovJumps(time, c);

        System.out.println("Reward conditional on X(0) = i, X(t) = j: " + new Vector(c));

        double endTime = 10.0;
        int steps = 10;
        for (time = 0.0; time < endTime; time += (endTime / steps)) {
            markovRewards.computeCondStatMarkovJumps(time, c);
            System.out.println(time + "," + c[0]);   // start = 0, end = 0
        }
    }

//    0 -> \alpha -> 1;   1 -> \beta -> 0
//    \eigenvectors =
//      \left[
//          1, -\alpha \\
//          1, \beta \\
//    \right]
//    \inveigenvectors =
//        \frac{1}{\alpha + \beta}
//     \left[
//            \beta,  \alpha, \\
//               -1, 1  \\
//        \right]

    private static double analyticProb(int from, int to, double alpha, double beta, double t) {

        double total = alpha + beta;
        int entry = from * 2 + to;

        switch (entry) {
            case 0: // 0 -> 0
                return (beta + alpha * Math.exp(-total * t)) / total;
            case 1: // 0 -> 1
                return (alpha - alpha * Math.exp(-total * t)) / total;
            case 2: // 1 -> 0
                return (beta - beta * Math.exp(-total * t)) / total;
            case 3: // 1 -> 1
                return (alpha + beta * Math.exp(-total * t)) / total;
            default:
                throw new RuntimeException();
        }
    }
}
