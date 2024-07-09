/*
 * SubstitutionRateMatrixMixtureParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.AminoAcidMixture;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.SubstitutionRateMatrixMixture;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class SubstitutionRateMatrixMixtureParser extends AbstractXMLObjectParser {
    private static final String MIXTURE_MODEL = "substitutionRateMatrixMixtureModel";

    public String getParserName() {
        return MIXTURE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<SubstitutionModel> modelList = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            SubstitutionModel model = (SubstitutionModel) xo.getChild(i);
            // TODO figure out why logic of this class does not apply to EAAM
            if (model instanceof EmpiricalAminoAcidModel) {throw new RuntimeException("For mixtures of empirical amino acid models use aminoAcidMixtureModel");}
            modelList.add(model);
        }

        return new SubstitutionRateMatrixMixture(modelList);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "List of Q-matrices for combining in a GLM.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
    };
}
