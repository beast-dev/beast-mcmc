/*
 * MG94CodonModel.java
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

package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.Codons;
import dr.inference.model.Parameter;

/**
 * Muse-Gaut model of codon evolution
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe lemey
 */
public class MG94CodonModel extends AbstractCodonModel {

    protected Parameter alphaParameter;
    protected Parameter betaParameter;

    private final int numSynTransitions;
    private final int numNonsynTransitions;

    public MG94CodonModel(Codons codonDataType, Parameter alphaParameter, Parameter betaParameter,
                          FrequencyModel freqModel) {
        this(codonDataType, alphaParameter, betaParameter, freqModel,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    public MG94CodonModel(Codons codonDataType,
                          Parameter alphaParameter,
                          Parameter betaParameter,
                          FrequencyModel freqModel, EigenSystem eigenSystem) {
        super("MG94", codonDataType, freqModel, eigenSystem);

        this.alphaParameter = alphaParameter;
        addVariable(alphaParameter);
        alphaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                alphaParameter.getDimension()));

        this.betaParameter = betaParameter;
        addVariable(betaParameter);
        betaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                betaParameter.getDimension()));

        numSynTransitions = getNumSynTransitions();
        numNonsynTransitions = getNumNonsynTransitions();
    }

    private int countRates(int i, int j) {
        int count = 0;

        for (byte rate : rateMap) {
            if (rate == i || rate == j) {
                count++;
            }
        }
        return count;
    }

    private int getNumSynTransitions() {
        return 2 * countRates(1, 2);
    }

    private int getNumNonsynTransitions() {
        return 2 * countRates(3, 4);
    }

    protected double getNormalizationValue(double[][] matrix, double[] pi) {
        double norm = 1.0;
        if (doNormalization) {
            double ratio =
                    getAlpha() + getBeta();
            //(numSynTransitions * getAlpha() + numNonsynTransitions * getBeta()) / (numSynTransitions + numNonsynTransitions);
            norm = super.getNormalizationValue(matrix, pi) / ratio;
        }
        return norm;
    }

    public double getAlpha() {
        return alphaParameter.getParameterValue(0);
    }

    public double getBeta() {
        return betaParameter.getParameterValue(0);
    }

    protected void setupRelativeRates(double[] rates) {

        double alpha = getAlpha() / numSynTransitions;
        double beta = getBeta() / numNonsynTransitions;
        for (int i = 0; i < rateCount; i++) {
            switch (rateMap[i]) {
                case 0:
                    rates[i] = 0.0;
                    break;            // codon changes in more than one codon position
                case 1:
                    rates[i] = alpha;
                    break;        // synonymous transition
                case 2:
                    rates[i] = alpha;
                    break;        // synonymous transversion
                case 3:
                    rates[i] = beta;
                    break;         // non-synonymous transition
                case 4:
                    rates[i] = beta;
                    break;            // non-synonymous transversion
            }
        }
    }

    public void setNormalization(boolean normalize) {
        this.doNormalization = normalize;
    }

    private boolean doNormalization = true;
}