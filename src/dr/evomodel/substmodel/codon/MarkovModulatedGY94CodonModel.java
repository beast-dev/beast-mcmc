/*
 * MarkovModulatedGY94CodonModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel.codon;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.HiddenCodons;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedGY94CodonModel extends GY94CodonModel {

    private static final byte RATE = 5;

    public MarkovModulatedGY94CodonModel(
            HiddenCodons codonDataType,
            Parameter switchingRates,
            Parameter omegaParameter,
            Parameter kappaParameter,
            FrequencyModel freqModel) {
        this(codonDataType, switchingRates, omegaParameter, kappaParameter, freqModel,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    public MarkovModulatedGY94CodonModel(
            HiddenCodons codonDataType,
            Parameter switchingRates,
            Parameter omegaParameter,
            Parameter kappaParameter,
            FrequencyModel freqModel,
            EigenSystem eigenSystem) {

        super(codonDataType, omegaParameter, kappaParameter, freqModel, eigenSystem);

        this.hiddenClassCount = codonDataType.getHiddenClassCount();
        this.switchingRates = switchingRates;
        addVariable(switchingRates);

        // Subclassed constructors fill relativeRates with 1
        for (int i = 0; i < relativeRates.length; i++)
            relativeRates[i] = 0.0;
    }

//    protected void handleVariableChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
//         // relativeRates changed
//        System.err.println("parameter "+parameter.getId()+" index = "+index+" changed");
//         updateMatrix = true;
//         ratesChanged();
////        fireModelChanged();
//     }


    protected void setupRelativeRates(double[] relativeRates) {
        double kappa = getKappa();
        double[] omega = omegaParameter.getParameterValues();
        double[] rates = switchingRates.getParameterValues();

        int stateCount = this.stateCount / hiddenClassCount;

        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                for (int h = 0; h < hiddenClassCount; h++) {

                    int d = getIndex(h * stateCount + i, h * stateCount + j, this.stateCount);
                    switch (rateMap[index]) {
                        case 0:
                            relativeRates[d] = 0.0;
                            break;                                // codon changes in more than one codon position
                        case 1:
                            relativeRates[d] = kappa;
                            break;                                // synonymous transition
                        case 2:
                            relativeRates[d] = 1.0;
                            break;                                // synonymous transversion
                        case 3:
                            relativeRates[d] = kappa * omega[h];
                            break;                                // non-synonymous transition
                        case 4:
                            relativeRates[d] = omega[h];
                            break;                                // non-synonymous transversion

                    }
                }
                index++;
            }
        }

        double[] freqs = freqModel.getFrequencies();
        int rateIndex = 0;
        for (int g = 0; g < hiddenClassCount; g++) {
            for (int h = g + 1; h < hiddenClassCount; h++) {  // from g -> h
                for (int i = 0; i < stateCount; i++) {
                    int d = getIndex(g * stateCount + i, h * stateCount + i, this.stateCount);
                    // correct for the fact that setupMatrix post-multiplies these rates
                    relativeRates[d] = rates[rateIndex] / freqs[i];
                }
                rateIndex++;
            }
        }
    }

    // Mapping: Matrix[i][j] = Compressed vector[i*(S - 3/2) - i^2 / 2 + j - 1]
    private static int getIndex(int i, int j, int S) {
        return (i * (2 * S - 3) - i * i) / 2 + j - 1;
    }

    protected void constructRateMap() {
        // Construct map for non-hidden states only
        hiddenClassCount = ((HiddenCodons) codonDataType).getHiddenClassCount();
        stateCount /= hiddenClassCount;
        super.constructRateMap();
        stateCount *= hiddenClassCount;
    }

    public static void main(String[] args) {
        GY94CodonModel codonModel = new GY94CodonModel(Codons.UNIVERSAL,
                new Parameter.Default(1.0), new Parameter.Default(2.0),
                new FrequencyModel(Codons.UNIVERSAL, new Parameter.Default(61, 1.0 / 61.0)));
        EigenDecomposition ed1 = codonModel.getEigenDecomposition();
//        double[][] q = codonModel.getQ();

//        System.err.println("matrixQ = \n"+codonModel.printQ());// new Matrix(q));
        FrequencyModel freqModel = new FrequencyModel(HiddenCodons.UNIVERSAL_HIDDEN_2, new Parameter.Default(122, 1.0 / 122.0));
        System.err.println("freq = " + new Vector(freqModel.getFrequencies()));
//        System.exit(-1);
        MarkovModulatedGY94CodonModel mmCodonModel = new MarkovModulatedGY94CodonModel(HiddenCodons.UNIVERSAL_HIDDEN_2,
                new Parameter.Default(2, 5.0), new Parameter.Default(2, 1.0),
                new Parameter.Default(2.0), freqModel
        );
        EigenDecomposition ed2 = mmCodonModel.getEigenDecomposition();
        System.err.println("matrixQ = \n" + mmCodonModel.printQ());// new Matrix(q));
    }

    protected double getMINFDIFF() {
        return 1.0E-10;
    }

    protected double getMINFREQ() {
        return 1.0E-10;
    }


    private int hiddenClassCount;
    private Parameter switchingRates;

}
