/*
 * AbstractCovarionDNAModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.OldHiddenNucleotides;
import dr.inference.model.Parameter;

/**
 * A general time reversible model of nucleotide evolution with
 * covarion hidden rate categories.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
abstract public class AbstractCovarionDNAModel extends AbstractSubstitutionModel {

    public static final String HIDDEN_CLASS_RATES = "hiddenClassRates";
    public static final String SWITCHING_RATES = "switchingRates";
    public static final String FREQUENCIES = "frequencies";

    /**
     * @param name             the name of the covarion substitution model
     * @param dataType         the data type
     * @param freqModel        the equlibrium frequencies
     * @param hiddenClassRates the relative rates of the hidden categories
     *                         (first hidden category has rate 1.0 so this parameter
     *                         has dimension one less than number of hidden categories.
     *                         each hidden category.
     * @param switchingRates   rate of switching between hidden categories
     */
    public AbstractCovarionDNAModel(String name,
                                    OldHiddenNucleotides dataType,
                                    Parameter hiddenClassRates,
                                    Parameter switchingRates,
                                    FrequencyModel freqModel) {

        super(name, dataType, freqModel);

        hiddenClassCount = dataType.getHiddenClassCount();

        this.hiddenClassRates = hiddenClassRates;
        this.switchingRates = switchingRates;

        assert hiddenClassRates.getDimension() == hiddenClassCount - 1;

        int hiddenClassCount = getHiddenClassCount();

        int switchingClassCount = hiddenClassCount * (hiddenClassCount - 1) / 2;

        if (switchingRates.getDimension() != switchingClassCount) {
            throw new IllegalArgumentException("switching rate parameter must have " +
                    switchingClassCount + " rates for " + hiddenClassCount + " classes");
        }
        addVariable(switchingRates);
        addVariable(hiddenClassRates);
        constructRateMatrixMap();
    }

    /**
     * @return the relative rates of A<->C, A<->G, A<->T, C<->G, C<->T and G<->T substitutions
     */
    abstract double[] getRelativeDNARates();

    /**
     * @return the number of hidden classes in this covarion model.
     */
    public final int getHiddenClassCount() {
        return hiddenClassCount;
    }

    public void frequenciesChanged() {
        // DO NOTHING
    }

    public void ratesChanged() {
        setupRelativeRates();
    }

    protected void setupRelativeRates() {

        double[] phi = switchingRates.getParameterValues();
        double[] rr = getRelativeDNARates();
        double[] hiddenRates = hiddenClassRates.getParameterValues();

        for (int i = 0; i < rateCount; i++) {
            if (rateMatrixMap[i] == 0) {
                relativeRates[i] = 0.0;
            } else if (rateMatrixMap[i] < 7) {
                if (hiddenClassMap[i] == 0) {
                    relativeRates[i] = rr[rateMatrixMap[i] - 1];
                } else {
                    relativeRates[i] = rr[rateMatrixMap[i] - 1] * hiddenRates[hiddenClassMap[i] - 1];
                }
            } else {
                relativeRates[i] = phi[rateMatrixMap[i] - 7];
            }
        }
    }

    /**
     * Construct a map of the rate classes in the rate matrix:
     * 0: class and nucleotide change
     * 1: A <-> C
     * 2: A <-> G
     * 3: A <-> T
     * 4: C <-> G
     * 5: C <-> T
     * 6: G <-> T
     * 7: 0 <-> 1 class change
     * 8: 0 <-> 2 class change
     * et cetera
     */
    private void constructRateMatrixMap() {

        byte rateClass;
        int fromNuc, toNuc;
        int fromRate, toRate;
        int count = 0;

        rateMatrixMap = new byte[rateCount];
        hiddenClassMap = new byte[rateCount];

        for (int i = 0; i < stateCount; i++) {

            for (int j = i + 1; j < stateCount; j++) {

                fromNuc = i % 4;
                toNuc = j % 4;
                fromRate = i / 4;
                toRate = j / 4;

                if (fromNuc == toNuc) {
                    // rate transition
                    if (fromRate == toRate) {
                        throw new RuntimeException("Shouldn't be possible");
                    }

                    rateClass = (byte) (7 + getIndex(fromRate, toRate, hiddenClassCount));
                } else if (fromRate != toRate) {
                    rateClass = 0;
                } else {
                    rateClass = (byte) (1 + getIndex(fromNuc, toNuc, 4));
                }

                rateMatrixMap[count] = rateClass;
                hiddenClassMap[count] = (byte) fromRate;
                count++;
            }
        }
    }

    private int getIndex(int from, int to, int size) {
        int index = 0;

        int f = from;
        while (f > 0) {
            index += size - 1;
            f -= 1;
            size -= 1;
        }
        index += to - from - 1;

        return index;
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     */
    void normalize(double[][] matrix, double[] pi) {
        double subst = 0.0;
        int dimension = pi.length;

        for (int i = 0; i < dimension; i++) {
            subst += -matrix[i][i] * pi[i];
        }

        // normalize, including switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst;
            }
        }

        double switchingProportion = 0.0;
        for (int i = 0; i < hiddenClassCount; i++) {
            for (int j = i + 1; j < hiddenClassCount; j++) {
                for (int l = 0; l < 4; l++) {

                    int x = i * 4 + l;
                    int y = j * 4 + l;

                    switchingProportion += matrix[x][y] * pi[y];
                    switchingProportion += matrix[y][x] * pi[x];
                }
            }
        }
        // normalize, removing switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / (1.0 - switchingProportion);
            }
        }
    }

    Parameter switchingRates;
    Parameter hiddenClassRates;

    byte[] rateMatrixMap;
    byte[] hiddenClassMap;

    private int hiddenClassCount;
}
