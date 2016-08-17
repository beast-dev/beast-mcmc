/*
 * CodonLabeling.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.Codons;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         An enum for different types of codon labelings in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */

public enum CodonLabeling {
    SYN("S"), // synonymous mutations
    NON_SYN("N"); // non-synonymous mutations

    CodonLabeling(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    public static double[] getRegisterMatrix(CodonLabeling labeling, Codons codonDataType) {
        return getRegisterMatrix(labeling, codonDataType, false);
    }

    public static double[] getRegisterMatrix(CodonLabeling labeling, Codons codonDataType, boolean base64) {
        final int stateCount = codonDataType.getStateCount();
        final int rateCount = ((stateCount - 1) * stateCount) / 2;

        byte[] rateMap = Codons.constructRateMap(
                rateCount,
                stateCount,
                codonDataType,
                codonDataType.getGeneticCode());

        double[] registerMatrix = new double[stateCount * stateCount];

        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                byte b = rateMap[index];
                if (
                        (labeling == SYN && (b == 1 || b == 2)) ||
                                (labeling == NON_SYN && (b == 3 || b == 4))
                        ) {
                    registerMatrix[j * stateCount + i] = registerMatrix[i * stateCount + j] = 1.0;
                }
                index++;
            }
        }

        if (base64) { // Expand matrix back out to the 4 x 4 x 4 stateSpace for a product chain
            double[] oldRegisterMatrix = registerMatrix;
            registerMatrix = new double[64 * 64];
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    if (oldRegisterMatrix[i * stateCount + j] == 1) {
                        registerMatrix[codonDataType.getCanonicalState(i) * 64
                                + codonDataType.getCanonicalState(j)] = 1;
                    }
                }
            }
        }
        return registerMatrix;
    }

    private final String text;

    public static CodonLabeling parseFromString(String text) {
        for (CodonLabeling scheme : CodonLabeling.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return null;
    }
}
