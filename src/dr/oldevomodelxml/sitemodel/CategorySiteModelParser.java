/*
 * CategorySiteModelParser.java
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

import dr.oldevomodel.sitemodel.CategorySiteModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class CategorySiteModelParser extends AbstractXMLObjectParser {

    public static final String SITE_MODEL = "categorySiteModel";
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String RATE_PARAMETER = "rates";
    public static final String CATEGORIES = "categories";
    public static final String CATEGORY_STATES = "states";
    public static final String CATEGORY_STRING = "values";
    public static final String RELATIVE_TO = "relativeTo";

    public String getParserName() {
        return SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(SUBSTITUTION_MODEL);
        SubstitutionModel substitutionModel = (SubstitutionModel) cxo.getChild(SubstitutionModel.class);

        cxo = xo.getChild(MUTATION_RATE);
        Parameter muParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RATE_PARAMETER);
        Parameter rateParam = null;
        int relativeTo = 0;
        if (cxo != null) {
            rateParam = (Parameter) cxo.getChild(Parameter.class);
            relativeTo = cxo.getIntegerAttribute(RELATIVE_TO);
        }

        cxo = xo.getChild(CATEGORIES);
        String categories = "";
        String states = "";
        if (cxo != null) {
            categories = cxo.getStringAttribute(CATEGORY_STRING);
            states = cxo.getStringAttribute(CATEGORY_STATES);
        }

        return new CategorySiteModel(substitutionModel, muParam, rateParam, categories, states, relativeTo);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a gamma distributed rates across sites";
    }

    public Class getReturnType() {
        return CategorySiteModel.class;
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
            new ElementRule(RATE_PARAMETER, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(RELATIVE_TO, true),
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(CATEGORIES, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(CATEGORY_STRING, true),
                    AttributeRule.newStringRule(CATEGORY_STATES, true)
            }, true),
    };

}
