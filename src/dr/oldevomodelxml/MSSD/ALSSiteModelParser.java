/*
 * ALSSiteModelParser.java
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

package dr.oldevomodelxml.MSSD;

import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * This site model is a copy of SiteModelParser explicitly so that it is not overridden
 * by the BEAGLE variant. This is because ALSTreeLikelihood requires a non-BEAGLE site
 * model.
 */
public class ALSSiteModelParser extends AbstractXMLObjectParser {

    public static final String ALS_SITE_MODEL = "alsSiteModel";

    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String SUBSTITUTION_RATE = "substitutionRate";
    public static final String RELATIVE_RATE = "relativeRate";
    public static final String GAMMA_SHAPE = "gammaShape";
    public static final String GAMMA_CATEGORIES = "gammaCategories";
    public static final String PROPORTION_INVARIANT = "proportionInvariant";


    public String[] getParserNames() {
        return new String[]{
                getParserName(), "beast_" + getParserName()
        };
    }

    public String getParserName() {
        return ALS_SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getElementFirstChild(SUBSTITUTION_MODEL);

        String msg = "";

        Parameter muParam = null;
        if (xo.hasChildNamed(SUBSTITUTION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else  if (xo.hasChildNamed(MUTATION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(MUTATION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else if (xo.hasChildNamed(RELATIVE_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(RELATIVE_RATE);

            msg += "\n  with initial relative rate = " + muParam.getParameterValue(0);
        }

        Parameter shapeParam = null;
        int catCount = 4;
        if (xo.hasChildNamed(GAMMA_SHAPE)) {
            final XMLObject cxo = xo.getChild(GAMMA_SHAPE);
            catCount = cxo.getIntegerAttribute(GAMMA_CATEGORIES);
            shapeParam = (Parameter) cxo.getChild(Parameter.class);

            msg += "\n  " + catCount + " category discrete gamma with initial shape = " + shapeParam.getParameterValue(0);
        }

        Parameter invarParam = null;
        if (xo.hasChildNamed(PROPORTION_INVARIANT)) {
            invarParam = (Parameter) xo.getElementFirstChild(PROPORTION_INVARIANT);
            msg += "\n  initial proportion of invariant sites = " + invarParam.getParameterValue(0);
        }

        Logger.getLogger("dr.evomodel").info("Creating site model." + (msg.length() > 0 ? msg : ""));

        return new GammaSiteModel(substitutionModel, muParam, shapeParam, catCount, invarParam);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a gamma distributed rates across sites";
    }

    public Class getReturnType() {
        return GammaSiteModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SubstitutionModel.class)
            }),
            new XORRule(
                    new XORRule(
                            new ElementRule(SUBSTITUTION_RATE, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
                    ),
                    new ElementRule(RELATIVE_RATE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }), true
            ),
            new ElementRule(GAMMA_SHAPE, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(GAMMA_CATEGORIES, true),
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(PROPORTION_INVARIANT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };
}
