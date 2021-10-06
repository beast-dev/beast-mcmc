/*
 * GLMSubstitutionModelParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.AminoAcidMixture;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class AminoAcidMixtureParser extends AbstractXMLObjectParser {

    public static final String GLM_SUBSTITUTION_MODEL = "aminoAcidMixtureModel";


    public String getParserName() {
        return GLM_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<EmpiricalAminoAcidModel> modeList = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            EmpiricalAminoAcidModel model = (EmpiricalAminoAcidModel) xo.getChild(i);
            modeList.add(model);
        }

        return new AminoAcidMixture(modeList);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Rates from named amino acid models for use by GLM substitution models.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(EmpiricalAminoAcidModel.class, 1, Integer.MAX_VALUE),
    };

}
