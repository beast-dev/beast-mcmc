/*
 * ProductChainFrequencyModelTest.java
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

package test.dr.evomodel.substmodel;

import test.dr.math.MathTestCase;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.ProductChainFrequencyModel;
import dr.evolution.datatype.Nucleotides;
import dr.math.matrixAlgebra.Vector;

import java.util.List;
import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */
public class ProductChainFrequencyModelTest extends MathTestCase {

    public void testFrequencyModel() {

        FrequencyModel firstPosition = new FrequencyModel(Nucleotides.INSTANCE, freq1);
        FrequencyModel secondPosition = new FrequencyModel(Nucleotides.INSTANCE, freq2);
        FrequencyModel thirdPosition = new FrequencyModel(Nucleotides.INSTANCE, freq3);

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>(3);
        freqModels.add(firstPosition);
        freqModels.add(secondPosition);
        freqModels.add(thirdPosition);

        ProductChainFrequencyModel pcFreqModel = new ProductChainFrequencyModel("freq", freqModels);

        double[] freqs = pcFreqModel.getFrequencies();

        System.out.println("Freq.length = " + freqs.length);

        int pos1 = 2;
        int pos2 = 1;
        int pos3 = 3;

        int index = computeIndex(pos1, pos2, pos3);
        
        System.out.println("Entry: " + new Vector(pcFreqModel.decomposeEntry(index)));
        System.out.println("Freq = " + freqs[index]);
        System.out.println("Freq = " + computeFreq(pos1, pos2, pos3));

        assertEquals(computeFreq(pos1, pos2, pos3), freqs[index]);       
    }

    private int computeIndex(int i, int j, int k) {
        return i * 16 + j * 4 + k;
    }

    private double computeFreq(int i, int j, int k) {
        return freq1[i] * freq2[j] * freq3[k];
    }

    private static double[] freq1 = {0.3, 0.25, 0.20, 0.25}; // A,C,G,T
    private static double[] freq2 = {0.1, 0.4,  0.25, 0.25}; // A,C,G,T
    private static double[] freq3 = {0.3, 0.15, 0.25, 0.30}; // A,C,G,T
}
