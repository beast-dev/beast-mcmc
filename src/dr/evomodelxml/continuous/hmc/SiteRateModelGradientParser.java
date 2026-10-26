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

import dr.evomodel.treedatalikelihood.PartitionedGradients;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.SiteRateModelGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses {@code <siteRateModelGradient>}: the analytic gradient of a tree data likelihood with respect to the gamma
 * shape, the proportion of invariant sites or the relative rate of its site model.
 * <p>
 * The likelihood may be a multi-partition {@code <treeDataLikelihood>} (a CompoundLikelihood of one tree data
 * likelihood per partition): a site-model parameter shared by several partitions gets the sum of their gradients,
 * and a compound of per-partition parameters, such as allNus, gets one block per partition. Give such a compound
 * its fixed sum with {@code simplexSum} (allNus sums to one), so gradient checks and transforms use the tangent.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SiteRateModelGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "siteRateModelGradient";
    private static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    private static final String SIMPLEX_SUM = "simplexSum";
    private static final double DEFAULT_TOLERANCE = 1E-4;

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final double tolerance = xo.getAttribute(GRADIENT_CHECK_TOLERANCE, DEFAULT_TOLERANCE);
        final double simplexSum = xo.getAttribute(SIMPLEX_SUM, Double.NaN);
        try {
            return PartitionedGradients.create(PartitionedGradients.treeDataLikelihoods(likelihood), parameter,
                    new PartitionedGradients.Factory() {
                        @Override
                        public boolean appliesTo(TreeDataLikelihood treeDataLikelihood, Parameter p) {
                            return SiteRateModelGradient.isSiteRateParameter(treeDataLikelihood, p);
                        }

                        @Override
                        public GradientWrtParameterProvider create(TreeDataLikelihood treeDataLikelihood, Parameter p) {
                            return new SiteRateModelGradient(treeDataLikelihood, p, tolerance);
                        }
                    }, simplexSum, tolerance);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(CompoundLikelihood.class)),
            new ElementRule(Parameter.class),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
            AttributeRule.newDoubleRule(SIMPLEX_SUM, true),
    };

    @Override
    public String getParserDescription() {
        return "Analytic gradient of a tree data likelihood (or a compound of per-partition ones) wrt the gamma shape,"
                + " proportion of invariant sites or relative rate of its site model(s)";
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }
}
