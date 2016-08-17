/*
 * MG94HKYCodonModel.java
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
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.Codons;
import dr.inference.model.Parameter;

/**
 * Muse-Gaut model of codon evolution
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe lemey
 */
public class MG94HKYCodonModel extends MG94CodonModel {

    protected Parameter kappaParameter;
//    protected Parameter AtoC_Parameter;
//    protected Parameter AtoG_Parameter;
//    protected Parameter AtoT_Parameter;
//    protected Parameter CtoG_Parameter;
//    protected Parameter CtoT_Parameter;
//    protected Parameter GtoT_Parameter;

    public MG94HKYCodonModel(Codons codonDataType, Parameter alphaParameter, Parameter betaParameter, Parameter kappaParameter,
                             FrequencyModel freqModel) {
        this(codonDataType, alphaParameter, betaParameter, kappaParameter, freqModel,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    public MG94HKYCodonModel(Codons codonDataType,
                             Parameter alphaParameter,
                             Parameter betaParameter,
                             Parameter kappaParameter,
//                          Parameter AtoC_Parameter,
//                          Parameter AtoG_Parameter,
//                          Parameter AtoT_Parameter,
//                          Parameter CtoG_Parameter,
//                          Parameter CtoT_Parameter,
//                          Parameter GtoT_Parameter,
                             FrequencyModel freqModel, EigenSystem eigenSystem) {
        super(codonDataType, alphaParameter, betaParameter, freqModel, eigenSystem);

        this.kappaParameter = kappaParameter;
        addVariable(kappaParameter);
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
        kappaParameter.getDimension()));

    }

    public double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    protected void setupRelativeRates(double[] rates) {

        double alpha = getAlpha() / numSynTransitions;
        double beta = getBeta() / numNonsynTransitions;

        double kappa = getKappa();
        for (int i = 0; i < rateCount; i++) {
            switch (rateMap[i]) {
                case 0:
                    rates[i] = 0.0;
                    break;            // codon changes in more than one codon position
                case 1:
                    rates[i] = alpha  * kappa;
                    break;        // synonymous transition
                case 2:
                    rates[i] = alpha;
                    break;        // synonymous transversion
                case 3:
                    rates[i] = beta  * kappa;
                    break;         // non-synonymous transition
                case 4:
                    rates[i] = beta;
                    break;            // non-synonymous transversion
            }
        }
    }

}