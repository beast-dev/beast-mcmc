/*
 * SiteRateModelGradientParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.SiteRateModelGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses {@code <siteRateModelGradient>}: the analytic gradient of a tree data likelihood with respect to the gamma
 * shape, the proportion of invariant sites or the relative rate of its site model.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SiteRateModelGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "siteRateModelGradient";
    private static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    private static final double DEFAULT_TOLERANCE = 1E-4;

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final double tolerance = xo.getAttribute(GRADIENT_CHECK_TOLERANCE, DEFAULT_TOLERANCE);
        try {
            return new SiteRateModelGradient(treeDataLikelihood, parameter, tolerance);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(Parameter.class),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
    };

    @Override
    public String getParserDescription() {
        return "Analytic gradient of a tree data likelihood wrt the gamma shape, proportion of invariant sites or"
                + " relative rate of its site model";
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }
}
