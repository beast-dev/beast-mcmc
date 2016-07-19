/*
 * SpecificEigendecompositionTest.java
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

import dr.evolution.datatype.GeneralDataType;
import dr.evomodel.substmodel.*;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

import java.util.Arrays;


public class SpecificEigendecompositionTest extends MathTestCase {

    public SpecificEigendecompositionTest() {
        super();
    }

    private ComplexSubstitutionModel setupModel(int dim, double[] rates) {
        String[] labels = new String[dim];
        for (Integer i = 0; i < dim; i++) {
            labels[i] = i.toString();
        }
        GeneralDataType dataType = new GeneralDataType(labels);
        Parameter freqVector = new Parameter.Default(dim);
        for (int i = 0; i < dim; i++) {
            freqVector.setParameterValue(i, (double) 1 / dim);
        }
        FrequencyModel freqModel = new FrequencyModel(dataType, freqVector);
        Parameter rateVector = new Parameter.Default(rates);

        return new ComplexSubstitutionModel("test", dataType, freqModel, rateVector) {
            protected EigenSystem getDefaultEigenSystem(int stateCount) {
                return new ComplexColtEigenSystem(stateCount, false, ColtEigenSystem.defaultMaxConditionNumber, ColtEigenSystem.defaultMaxIterations);
            }
        };
    }

//    private static int dim = 4;
//
//    private static double[] testRates = {
//     0.27577,     0.39669,     0.07341,
//                  0.07491,     0.08690,
//                                   0.0,
//     0.75, 0.0, 0.19029,
//           0.0, 0.44924,
//                0.70278
//    };
//
//    private static double[] checkEigenvalues = {-1.7950052, -1.7884673, -0.4165275,  0.0000000, 0.0, 0.0, 0.0, 0.0};

//    private static int dim = 5;
//
//    private static double[] testRates = {
//            0.0,                      0.0,          0.0,                    0.0,
//                                      0.0,          0.0,                    0.0,
//                                                    3.478773070125323,      0.2848265288341367,
//                                                                            0.0,
//    0.0,    3.606696517465288,        0.0,          0.3708731237136557,
//            4.421855313152474,        0.0,          2.6121833491266533,
//                                      0.0,          3.997241604528838,
//                                                    1.227550493053631
//    };
//
//    private static double[] checkEigenvalues = {-3.0214359488392066, -1.978564049686774, 0.0, 0.0, 0.0};

    private static int dim = 6;

    private static double[] testRates = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0,
            3.5409959027854936,
            0.0, 0.0, 0.0, 1.9115985677138407, 3.6880365950035827,
            0.0, 0.0, 3.135364040895322, 1.2936152589276488,
            0.0, 2.4748026323565226, 4.738559654533422,
            3.4092727663178453, 2.1723540463137088,
            3.635400535152609

    };

    private static double[] checkEigenvalues = {-3.72531, -2.27469, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


//    private static double[] testRates = {
//        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
//    };


    private static double tolerance = 1E-4;

    public void testEigendecomposition() {
        System.out.println("Testing specific eigendecomposition...");
        ComplexSubstitutionModel csm = setupModel(dim, testRates);
        double[] tmp = new double[dim * dim];
        csm.getInfinitesimalMatrix(tmp);
        System.out.println("Rates: " + new Vector(tmp) + "\n");

        EigenDecomposition e = csm.getEigenDecomposition();
        System.out.println("Val: " + new Vector(e.getEigenValues()));
        System.out.println("Vec: " + new Vector(e.getEigenVectors()));
        System.out.println("Inv: " + new Vector(e.getInverseEigenVectors()));

        csm.getTransitionProbabilities(1.0, tmp);
        System.out.println(new Vector(tmp));

        double[] eigenValues = e.getEigenValues();
        Arrays.sort(eigenValues);

        assertEquals(checkEigenvalues, e.getEigenValues(), tolerance);

    }

}
