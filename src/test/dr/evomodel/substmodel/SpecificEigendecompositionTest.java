/*
 * SpecificEigendecompositionTest.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beagle.evomodel.substmodel.ComplexSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.GeneralDataType;
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

        return new ComplexSubstitutionModel("test", dataType, freqModel, rateVector);
    }

    private static double[] testRates = {
     0.27577,     0.39669,     0.07341,
                  0.07491,     0.08690,
                                   0.0,
     0.75, 0.0, 0.19029,
           0.0, 0.44924,
                0.70278
    };

//    private static double[] testRates = {
//        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
//    };


    private static double[] checkEigenvalues = {-1.7950052, -1.7884673, -0.4165275,  0.0000000, 0.0, 0.0, 0.0, 0.0};

    private static double tolerance = 1E-4;

    public void testEigendecomposition() {
        System.out.println("Testing specific eigendecomposition...");
        ComplexSubstitutionModel csm = setupModel(4, testRates);
        double[] tmp = new double[16];
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
