/*
 * SiteRateCategoryGibbsOperatorParser.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodelxml.treedatalikelihood.MAT;

import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.MATDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.SiteRateCategoryGibbsOperator;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Parser for the Gibbs operator on the site rate categories of a MAT likelihood.
 */
public class SiteRateCategoryGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String SITE_RATE_CATEGORY_GIBBS_OPERATOR = "siteRateCategoryGibbsOperator";

    @Override
    public String getParserName() {
        return SITE_RATE_CATEGORY_GIBBS_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        DataLikelihoodDelegate delegate = likelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof MATDataLikelihoodDelegate)) {
            throw new XMLParseException("The likelihood of " + xo.getId() + " must be a MAT likelihood (matDataLikelihood)");
        }
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        try {
            return new SiteRateCategoryGibbsOperator((MATDataLikelihoodDelegate) delegate, weight);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage() + " in " + xo.getId());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Gibbs update of the site rate categories of a MAT likelihood whose site model has more than one " +
                "rate category: each variable site's category is drawn from its full conditional.";
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Class getReturnType() {
        return SiteRateCategoryGibbsOperator.class;
    }
}
