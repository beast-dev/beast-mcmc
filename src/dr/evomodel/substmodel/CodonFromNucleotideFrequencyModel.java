/*
 * CodonFromNucleotideFrequencyModel.java
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
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class CodonFromNucleotideFrequencyModel extends FrequencyModel {

    private final FrequencyModel nucleotideFrequencyModel;
    private final Codons dataType;

    public CodonFromNucleotideFrequencyModel(Codons dataType,
                                             FrequencyModel nucleotideFrequencyModel,
                                             Parameter codonFrequencies) {
        super(dataType, codonFrequencies);

        if (nucleotideFrequencyModel.getDataType() != Nucleotides.INSTANCE) {
            throw new IllegalArgumentException("Must provide a nucleotide frequency model");
        }

        this.dataType = dataType;
        this.nucleotideFrequencyModel = nucleotideFrequencyModel;
        updateFrequencyParameter();
        addModel(nucleotideFrequencyModel);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        assert (model == nucleotideFrequencyModel);

        updateFrequencyParameter();

        fireModelChanged(model);
    }

//    public void setFrequency(int i, double value) {
//        throw new RuntimeException("Should not directly set frequencies");
//    }

    private void updateFrequencyParameter() {

        // TODO Maybe make this lazy if necessary

        for (int nuc1 = 0; nuc1 < 4; ++nuc1) {

            double freq1 = nucleotideFrequencyModel.getFrequency(nuc1);

            for (int nuc2 = 0; nuc2 < 4; ++nuc2) {

                double freq2 = nucleotideFrequencyModel.getFrequency(nuc2);

                for (int nuc3 = 0; nuc3 < 4; ++nuc3) {

                    double freq3 = nucleotideFrequencyModel.getFrequency(nuc3);

                    final int state = dataType.getState(nuc1, nuc2, nuc3);

                    if (!(dataType.isStopCodon(state))) {

                        frequencyParameter.setParameterValue( state,
                                freq1 * freq2 * freq3);

                    }
                }
            }
        }

        final double sum = getSumOfFrequencies(frequencyParameter);

        for (int i = 0; i < frequencyParameter.getSize(); i++) {

            frequencyParameter.setParameterValue(i, frequencyParameter.getParameterValue(i) / sum);

        }
    }
}
