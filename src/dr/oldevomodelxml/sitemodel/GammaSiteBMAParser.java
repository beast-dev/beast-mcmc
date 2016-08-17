/*
 * GammaSiteBMAParser.java
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

package dr.oldevomodelxml.sitemodel;

import dr.xml.*;
import dr.oldevomodel.sitemodel.GammaSiteBMA;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for Gamma site BMA.
 *
 */
public class GammaSiteBMAParser extends AbstractXMLObjectParser {

    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String LOG_GAMMA_SHAPE = "logGammaShape";
    public static final String GAMMA_CATEGORIES = "gammaCategories";
    public static final String LOGIT_PROPORTION_INVARIANT = "logitProportionInvariant";
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String GAMMA_SITE_BMA = "gammaSiteBMA";

    public String getParserName() {
        return GAMMA_SITE_BMA;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getElementFirstChild(SUBSTITUTION_MODEL);

        Parameter muParam = (Parameter) xo.getElementFirstChild(MUTATION_RATE);

        Parameter logitInvar = (Parameter) xo.getElementFirstChild(LOGIT_PROPORTION_INVARIANT);

        final XMLObject cxo = xo.getChild(LOG_GAMMA_SHAPE);
        Parameter logShape = (Parameter) cxo.getChild(Parameter.class);

        int catCount = cxo.getIntegerAttribute(GAMMA_CATEGORIES);

        Variable<Integer> modelChoose = (Variable<Integer>) xo.getElementFirstChild(MODEL_CHOOSE);

        return new GammaSiteBMA(substitutionModel, muParam, logitInvar, logShape, catCount, modelChoose);
    }

    public String getParserDescription() {
        return "A SiteModel that does BMA for a gamma distributed rates across sites";
    }

    public Class getReturnType() {
        return GammaSiteBMA.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SubstitutionModel.class)
            }),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(LOG_GAMMA_SHAPE, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(GAMMA_CATEGORIES),
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(LOGIT_PROPORTION_INVARIANT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(MODEL_CHOOSE, new XMLSyntaxRule[]{
                    new ElementRule(Variable.class)
            })
    };

}
