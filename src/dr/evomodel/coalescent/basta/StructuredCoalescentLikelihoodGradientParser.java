/*
 * StructuredCoalescentLikelihoodGradientParser.java
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

package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

public class StructuredCoalescentLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "structuredCoalescentLikelihoodGradient";
    private static final String WRT_PARAMETER = "wrtParameter";

//    @Override
//    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//        StructuredCoalescentLikelihood structuredCoalescentLikelihood = (StructuredCoalescentLikelihood) xo.getChild(StructuredCoalescentLikelihood.class);
//        GlmSubstitutionModel substitutionModel = (GlmSubstitutionModel) xo.getChild(GlmSubstitutionModel.class);
//        String wrtParamter = (String) xo.getAttribute(WRT_PARAMETER);
//
//        StructuredCoalescentLikelihoodGradient.WrtParameter type = StructuredCoalescentLikelihoodGradient.WrtParameter.factory(wrtParamter);
//
//        return new StructuredCoalescentLikelihoodGradient(structuredCoalescentLikelihood, substitutionModel, type);
//    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        BastaLikelihood BastaLikelihood = (BastaLikelihood) xo.getChild(BastaLikelihood.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        String wrtParamter = (String) xo.getAttribute(WRT_PARAMETER);

        StructuredCoalescentLikelihoodGradient.WrtParameter type = StructuredCoalescentLikelihoodGradient.WrtParameter.factory(wrtParamter);

        return new StructuredCoalescentLikelihoodGradient(BastaLikelihood, substitutionModel, type);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(BastaLikelihood.class),
            new ElementRule(SubstitutionModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return StructuredCoalescentLikelihoodGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
