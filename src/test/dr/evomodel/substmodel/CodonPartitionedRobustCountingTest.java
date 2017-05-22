/*
 * CodonPartitionedRobustCountingTest.java
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

import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.*;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class CodonPartitionedRobustCountingTest extends TestCase {

    private Codons codons;
    private GeneticCode geneticCode;
    private ProductChainSubstitutionModel productChainModel;

    public void setUp() {
        codons = Codons.UNIVERSAL;
        geneticCode = codons.getGeneticCode();
    }

    private void setUpNucleotides() {

        FrequencyModel freqModel0 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});
        FrequencyModel freqModel1 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});
        FrequencyModel freqModel2 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});

        HKY baseModel0 = new HKY(2.0, freqModel0);
        HKY baseModel1 = new HKY(2.0, freqModel1);
        HKY baseModel2 = new HKY(2.0, freqModel2);

        List<SubstitutionModel> baseModels = new ArrayList<SubstitutionModel>(3);
        baseModels.add(baseModel0);
        baseModels.add(baseModel1);
        baseModels.add(baseModel2);

        productChainModel = new ProductChainSubstitutionModel("productChain", baseModels);

        /*
            hky0 = as.eigen(hky.model(2, 1, c(0.25,0.25, 0.25, 0.25), scale = T))
            hky1 = as.eigen(hky.model(2, 1, c(0.25,0.25, 0.25, 0.25), scale = T))
            hky2 = as.eigen(hky.model(2, 1, c(0.25,0.25, 0.25, 0.25), scale = T))
            pc = ind.codon.eigen(hky0, hky1, hky2)
            syn.matrix = regist.synonym()

         */
    }

    private double sumMatrix(double[] mat) {
        double total = 0;
        for (double x : mat) {
            total += x;
        }
        return total;
    }

    public void testSumRegisterMatrices() {

        System.out.println("Testing the sum of the register matrices...");

        double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, true);
        double totalSyn = sumMatrix(synRegMatrix);
        System.out.println("Sum syn matrix = " + totalSyn);
        assertEquals(134, (int) totalSyn);

        /* R markovjumps
        sum(regist.synonym())
         */

        double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, true);
        double totalNonSyn = sumMatrix(nonSynRegMatrix);
        System.out.println("Sum non-syn matrix = " + totalNonSyn);
        assertEquals(392, (int) totalNonSyn);

        /* R markovjumps
        sum(regist.nonsynonym())
         */

        System.out.println("");
    }

    public void testStopCodons() {
        // AGA -> TGA
        // AAA -> TAA

        System.out.println("Testing several stop codons...");

        setUpNucleotides();

        int from, to;
        double[] synRegMatrix, nonSynRegMatrix;

        from = codons.getState('A', 'G', 'A');
        to = codons.getState('T', 'G', 'A');

        synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, true);
        nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, true);

        System.out.println("AGA -> TGA syn = " + synRegMatrix[from * 64 + to]);
        System.out.println("AGA -> TGA non = " + nonSynRegMatrix[from * 64 + to]);

        from = codons.getState('A', 'A', 'A');
        to = codons.getState('T', 'A', 'A');

        synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, true);
        nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, true);

        System.out.println("AAA -> TAA syn = " + synRegMatrix[from * 64 + to]);
        System.out.println("AAA -> TAA non = " + nonSynRegMatrix[from * 64 + to]);
    }

    public void testCounts() {

        System.out.println("Testing several conditional means...");

        setUpNucleotides();

        int from = codons.getState('A', 'A', 'A');
        int to = codons.getState('A', 'A', 'C');
        double time = 1.0;

        // Syn check
        MarkovJumpsSubstitutionModel markovJumps = new MarkovJumpsSubstitutionModel(productChainModel);
        double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, true);
        markovJumps.setRegistration(synRegMatrix);

        double[] c = new double[64 * 64];
        markovJumps.computeCondStatMarkovJumps(time, c);

        System.out.println("Cond 1 = " + c[from * 64 + to]);
        assertEquals(0.525401, c[from * 64 + to], 1E-4);

        to = codons.getState('A', 'A', 'G');
        System.out.println("Cond 2 = " + c[from * 64 + to]);
        assertEquals(0.970012, c[from * 64 + to], 1E-4);

        // Nonsyn check
        double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, true);
        markovJumps.setRegistration(nonSynRegMatrix);
        markovJumps.computeCondStatMarkovJumps(time, c);

        System.out.println("Cond 3 = " + c[from * 64 + to]);
        assertEquals(1.099807, c[from * 64 + to], 1E-4);

        System.out.println("");

        /*
            cm.syn    = cond.mean.markov.jumps(pc, regist.synonym(),    1.0)
            cm.nonsyn = cond.mean.markov.jumps(pc, regist.nonsynonym(), 1.0)
         */
    }

    public void testRegistrationMatrixCodonBase() {
        runRegistrationMatrixTest(false); // test in 64 - stop codon states
    }

    public void testRegistrationMatrix64Base() {
        runRegistrationMatrixTest(true); // test in 64 states (product chain)
    }

    public void runRegistrationMatrixTest(boolean base64) {

        System.out.println("Testing registration matrices...");

        byte[] rateMat = Codons.constructRateMap(codons);
        int stateCount = codons.getStateCount();

        if (base64) {
            rateMat = expandToBase64(codons, rateMat);
            stateCount = 64;
        }

        double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, base64);
        double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, base64);

        // Check all pair-wise code states
        int index = 0;
        boolean passed = true;
        for (int i = 0; i < stateCount; i++) {

            int stateI;
            if (base64) {
                stateI = i;
            } else {
                stateI = codons.getCanonicalState(i);
            }
            char iAA = geneticCode.getAminoAcidChar(stateI);

            for (int j = i + 1; j < stateCount; j++) {
                if (rateMat[index] != 0) {

                    int stateJ;
                    if (base64) {
                        stateJ = j;
                    } else {
                        stateJ = codons.getCanonicalState(j);
                    }
                    char jAA = geneticCode.getAminoAcidChar(stateJ);

                    if (iAA == jAA) { // Syn
                        if (!compare(synRegMatrix, nonSynRegMatrix, i, j, iAA, jAA, 1, 0, stateCount)) {
                            passed = false;
                        }
                    } else { // Non-syn
                        if (!compare(synRegMatrix, nonSynRegMatrix, i, j, iAA, jAA, 0, 1, stateCount)) {
                            passed = false;
                        }
                    }
                } else { // 0 rate
                    if (synRegMatrix[i * stateCount + j] != 0 || nonSynRegMatrix[i * stateCount + j] != 0 ||
                            synRegMatrix[j * stateCount + i] != 0 || nonSynRegMatrix[j * stateCount + i] != 0) {
                        System.out.println("fail on 0-rate " + i + " " + j);
                        passed = false;
                    }
                }
                index++;
            }
        }
        System.out.println("Passed " + stateCount + ": " + (passed ? "true" : "false") + "\n");
        assertTrue(passed);
    }

    private boolean compare(double[] synRegMatrix, double[] nonSynRegMatrix,
                            int i, int j,
                            char iAA, char jAA,
                            int value0, int value1,
                            int stateCount) {

        boolean passed = true;
        if ((synRegMatrix[i * stateCount + j] == value0 && nonSynRegMatrix[i * stateCount + j] == value1) &&

                (synRegMatrix[j * stateCount + i] == value0 && nonSynRegMatrix[j * stateCount + i] == value1)) {
//            System.out.print("compare: " + iAA + jAA + " ");
//            System.out.println("pass");
        } else {
            System.out.print("compare: " + iAA + jAA + " ");
            System.out.println("fail " + i + " " + j + " " + synRegMatrix[i * stateCount + j]
                    + " " + nonSynRegMatrix[i * stateCount + j]);
            passed = false;
        }
        return passed;
    }

    private byte[] expandToBase64(Codons codonDataType, byte[] in) {
        final int stateCount = codonDataType.getStateCount();
        byte[] intermediateOut = new byte[64 * 64];
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                intermediateOut[codonDataType.getCanonicalState(i) * 64
                        + codonDataType.getCanonicalState(j)] = in[index];
                index++;
            }
        }

        byte[] out = new byte[(64 * 63) / 2];
        index = 0;
        for (int i = 0; i < 64; i++) {
            for (int j = i + 1; j < 64; j++) {
                out[index] = intermediateOut[i * 64 + j];
                index++;
            }
        }
        return out;
    }
}
