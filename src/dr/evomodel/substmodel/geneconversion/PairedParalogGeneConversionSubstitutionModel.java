/*
 * TwoParalogGeneConversionSubstitutionModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel.geneconversion;

import dr.evolution.datatype.PairedDataType;
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogGeneConversionSubstitutionModel extends BaseSubstitutionModel implements Citable {

    private final BaseSubstitutionModel baseSubstitutionModel;
    private final Parameter igcRateParameter;

    public PairedParalogGeneConversionSubstitutionModel(String name,
                                                        BaseSubstitutionModel baseSubstitutionModel,
                                                        Parameter relativeGeneConversionRateParameter,
                                                        PairedDataType dataType) {

        super(name, dataType, new PairedParalogFrequencyModel(dataType, baseSubstitutionModel.getFrequencyModel().getFrequencyParameter()));
        this.baseSubstitutionModel = baseSubstitutionModel;
        this.igcRateParameter = relativeGeneConversionRateParameter;


    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using igc extension model for gene conversion rate estimations.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author( "A", "Griffing"),
                    new Author("J", "Thorne"),
            },
            "A phylogenetic approach finds abundant interlocus gene conversion in yeast",
            "Molecular Biology and Evolution",
            Citation.Status.PUBLISHED);

    @Override
    protected void frequenciesChanged() {

    }

    @Override
    protected void ratesChanged() {

    }

    @Override
    protected void setupRelativeRates(double[] rates) {

    }
}
