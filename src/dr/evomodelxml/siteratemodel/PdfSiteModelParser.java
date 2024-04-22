/*
 * PdfSiteModelParser.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.siteratemodel;

import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.PdfSiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class PdfSiteModelParser extends AbstractXMLObjectParser {

    private static final String SITE_MODEL = "pdfSiteModel";
    private static final String SUBSTITUTION_MODEL = "substitutionModel";
    private static final String RATES = "rates";
    private static final String WEIGHTS = "weights";

    public String getParserName() {
        return SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rates = (Parameter) xo.getElementFirstChild(RATES);
        Parameter weights = (Parameter) xo.getElementFirstChild(WEIGHTS);

        PdfSiteRateModel siteRateModel = new PdfSiteRateModel(SITE_MODEL, rates, weights);

        if (xo.hasChildNamed(SUBSTITUTION_MODEL)) {
            siteRateModel.setSubstitutionModel(
                    (SubstitutionModel) xo.getElementFirstChild(SUBSTITUTION_MODEL));

        }

        return siteRateModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteRateModel that has probability free distributed rates across sites";
    }

    public Class<GammaSiteRateModel> getReturnType() {
        return GammaSiteRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(WEIGHTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SubstitutionModel.class)
            }, true),
    };
}
