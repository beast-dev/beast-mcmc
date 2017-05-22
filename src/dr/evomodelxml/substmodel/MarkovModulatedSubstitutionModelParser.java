/*
 * MarkovModulatedSubstitutionModelParser.java
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

import dr.evolution.datatype.HiddenDataType;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.MarkovModulatedSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String MARKOV_MODULATED_MODEL = "markovModulatedSubstitutionModel";
    //    public static final String HIDDEN_COUNT = "hiddenCount";
    public static final String SWITCHING_RATES = "switchingRates";
    //    public static final String DIAGONALIZATION = "diagonalization";
    public static final String RATE_SCALAR = "rateScalar";
    public static final String GEOMETRIC_RATES = "geometricRates";
    public static final String RENORMALIZE = "renormalize";

    public String getParserName() {
        return MARKOV_MODULATED_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);
        if (!(dataType instanceof HiddenDataType)) {
            throw new XMLParseException("Must construct " + MARKOV_MODULATED_MODEL + " with hidden data types");
        }

        Parameter switchingRates = (Parameter) xo.getElementFirstChild(SWITCHING_RATES);

        List<SubstitutionModel> substModels = new ArrayList<SubstitutionModel>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object cxo = xo.getChild(i);
            if (cxo instanceof SubstitutionModel) {
                substModels.add((SubstitutionModel) cxo);
            }
        }

        boolean geometricRates = xo.getAttribute(GEOMETRIC_RATES, false);

        Parameter rateScalar = xo.hasChildNamed(RATE_SCALAR) ?
                (Parameter) xo.getChild(RATE_SCALAR).getChild(Parameter.class) : null;

        SiteRateModel siteRateModel = (SiteRateModel) xo.getChild(SiteRateModel.class);
        if (siteRateModel != null) {
            if (siteRateModel.getCategoryCount() != substModels.size() &&  substModels.size() % siteRateModel.getCategoryCount()  != 0) {
                throw new XMLParseException(
                        "Number of gamma categories must equal number of substitution models in " + xo.getId());
            }
        }

        MarkovModulatedSubstitutionModel mmsm = new MarkovModulatedSubstitutionModel(xo.getId(), substModels, switchingRates, dataType, null,
                rateScalar, geometricRates, siteRateModel);

        if (xo.getAttribute(RENORMALIZE, false)) {
            mmsm.setNormalization(true);
        }

        return mmsm;
    }

    public String getParserDescription() {
        return "This element represents the a Markov-modulated substitution model.";
    }

    public Class getReturnType() {
        return MarkovModulatedSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            AttributeRule.newStringRule(DataType.DATA_TYPE),
//            AttributeRule.newStringRule(GeneticCode.GENETIC_CODE),
//            new ElementRule(OMEGA,
//                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//            new ElementRule(KAPPA,
//                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SWITCHING_RATES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//            new ElementRule(FrequencyModel.class),
//            AttributeRule.newStringRule(DIAGONALIZATION),
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(GEOMETRIC_RATES, true),
            AttributeRule.newBooleanRule(RENORMALIZE, true),
            new ElementRule(RATE_SCALAR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),

            new ElementRule(SiteRateModel.class, true),
    };
}
