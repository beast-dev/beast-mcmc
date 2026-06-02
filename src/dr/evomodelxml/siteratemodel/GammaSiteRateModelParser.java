/*
 * GammaSiteRateModelParser.java
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

package dr.evomodelxml.siteratemodel;

import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateDelegate;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * This is a replacement to GammaSiteModelParser that uses the modular
 * DiscretizedSiteRateModel with a Gamma delegate.
 * @author Andrew Rambaut
 */
public class GammaSiteRateModelParser extends AbstractXMLObjectParser {

    public static final String GAMMA_SITE_RATE_MODEL = "gammaSiteRateModel";
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String SUBSTITUTION_RATE = "substitutionRate";
    public static final String RELATIVE_RATE = "relativeRate";
    public static final String WEIGHT = "weight";
    public static final String GAMMA_SHAPE = "gammaShape";
    public static final String CATEGORIES = "categories";
    public static final String PROPORTION_INVARIANT = "proportionInvariant";

    public String getParserName() {
        return GAMMA_SITE_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String msg = "";
        SubstitutionModel substitutionModel = null;

        double muWeight = 1.0;

        Parameter muParam = null;
        if (xo.hasChildNamed(SUBSTITUTION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else  if (xo.hasChildNamed(MUTATION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(MUTATION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else if (xo.hasChildNamed(RELATIVE_RATE)) {
            XMLObject cxo = xo.getChild(RELATIVE_RATE);
            muParam = (Parameter) cxo.getChild(Parameter.class);
            msg += "\n  with initial relative rate = " + muParam.getParameterValue(0);
            if (cxo.hasAttribute(WEIGHT)) {
                muWeight = cxo.getDoubleAttribute(WEIGHT);
                msg += " with weight: " + muWeight;
            }
        }
        
        if(xo.hasChildNamed(SUBSTITUTION_MODEL)){
            substitutionModel = (SubstitutionModel)xo.getElementFirstChild(SUBSTITUTION_MODEL);
        }

        int catCount = 4;
        catCount = xo.getIntegerAttribute(CATEGORIES);

        Parameter shapeParam = null;
        if (xo.hasChildNamed(GAMMA_SHAPE)) {
            XMLObject cxo = xo.getChild(GAMMA_SHAPE);

            shapeParam = (Parameter) cxo.getChild(Parameter.class);

            msg += "\n  " + catCount + " category discrete gamma with initial shape = " + shapeParam.getParameterValue(0);
            msg += "\n  using equal weight discretization of gamma distribution";
        }

        Parameter invarParam = null;
        if (xo.hasChildNamed(PROPORTION_INVARIANT)) {
            invarParam = (Parameter) xo.getElementFirstChild(PROPORTION_INVARIANT);
            msg += "\n  initial proportion of invariant sites = " + invarParam.getParameterValue(0);
        }

        if (msg.length() > 0) {
            Logger.getLogger("dr.evomodel").info("\nCreating site rate model: " + msg);
        } else {
            Logger.getLogger("dr.evomodel").info("\nCreating site rate model.");
        }

        GammaSiteRateDelegate delegate = new GammaSiteRateDelegate("GammaSiteRateDelegate", shapeParam, catCount, invarParam);

        DiscretizedSiteRateModel siteRateModel = new DiscretizedSiteRateModel(SiteModel.SITE_MODEL, muParam, muWeight, delegate);

        if (substitutionModel != null) {
            siteRateModel.setSubstitutionModel(substitutionModel);
            siteRateModel.addModel(substitutionModel);
        }
        
        return siteRateModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A DiscretizedSiteRateModel that has a gamma distributed rates across sites";
    }

    @Override
    public String[] getParserNames() {
        return super.getParserNames();
    }

    public Class<DiscretizedSiteRateModel> getReturnType() {
        return DiscretizedSiteRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

            AttributeRule.newIntegerRule(CATEGORIES, false),
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
                            AttributeRule.newDoubleRule(WEIGHT, true),
                            new ElementRule(Parameter.class)
                    }), true
            ),

            new ElementRule(GAMMA_SHAPE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),

            new ElementRule(PROPORTION_INVARIANT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };

}//END: class
