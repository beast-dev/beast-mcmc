/*
 * SquaredSplineCoalescentLikelihoodParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.coalescent;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.SquaredSplineCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * XML parser for SquaredSplineCoalescentLikelihood.
 *
 * The <rateParameter> element must be a CompoundParameter=[coefficients, intercept] defined
 * in XML and shared with the HMC operator and any prior gradient wrapper.  Requiring it
 * explicitly prevents the anonymous-internal-parameter anti-pattern: if the likelihood
 * created the compound parameter internally, the object returned by getParameter() would
 * be invisible to XML, and combining it with a prior gradient via JointGradient would fail
 * because JointGradient checks parameter-object equality across providers.
 *
 * Example XML:
 *
 * <compoundParameter id="squaredSplineRate">
 *     <parameter idref="spline.coefficients"/>
 *     <parameter idref="spline.intercept"/>
 * </compoundParameter>
 *
 * <squaredSplineCoalescentLikelihood>
 *
 *     <!-- one or more trees for multilocus analysis -->
 *     <treeModel idref="treeModel1"/>
 *
 *     <splines>
 *         <IntegratedTransformedSplines idref="splineRates"/>
 *     </splines>
 *
 *     <!-- the rate parameter used by HMC and gradient providers -->
 *     <rateParameter>
 *         <compoundParameter idref="squaredSplineRate"/>
 *     </rateParameter>
 *
 *     <epsilon>
 *         <parameter id="epsilon" value="1e-6" lower="0.0"/>
 *     </epsilon>
 *
 *     <!-- optional; default 1.0 per tree -->
 *     <ploidy>
 *         <parameter value="1.0 1.0"/>
 *     </ploidy>
 *
 * </squaredSplineCoalescentLikelihood>
 */
public class SquaredSplineCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME   = "squaredSplineCoalescentLikelihood";
    public static final String SPLINES       = "splines";
    public static final String RATE_PARAM    = "rateParameter";
    public static final String EPSILON       = "epsilon";
    public static final String PLOIDY        = "ploidy";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // --- splines ---
        XMLObject splinesXO = xo.getChild(SPLINES);
        SquaredCachedSplines splines =
                (SquaredCachedSplines) splinesXO.getChild(SquaredCachedSplines.class);

        // --- rateParameter: the shared compound [coefficients, intercept] ---
        // Must be defined in XML and passed here so the same object is used by
        // the likelihood, the HMC operator, and any prior gradient wrapper.
        XMLObject rateParamXO = xo.getChild(RATE_PARAM);
        Parameter rateParameter = (Parameter) rateParamXO.getChild(Parameter.class);

        // Validate ordering matches getGradientWrtParameters: [coefficients..., intercept]
        Parameter coefficients = splines.getCoefficients();
        Parameter intercept    = splines.getIntercept();
        if (rateParameter instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) rateParameter;
            if (cp.getDimension() != coefficients.getDimension() + intercept.getDimension()) {
                throw new XMLParseException(PARSER_NAME + " rateParameter dimension mismatch");
            }
        }

        // --- trees: accept TreeModel children or pre-built BigFastTreeIntervals ---
        List<BigFastTreeIntervals> intervalsList = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof TreeModel) {
                intervalsList.add(new BigFastTreeIntervals((TreeModel) child));
            } else if (child instanceof BigFastTreeIntervals) {
                intervalsList.add((BigFastTreeIntervals) child);
            }
        }

        if (intervalsList.isEmpty()) {
            throw new XMLParseException(
                    PARSER_NAME + ": at least one treeModel or bigFastTreeIntervals is required");
        }

        // --- epsilon ---
        Parameter epsilon;
        if (xo.hasChildNamed(EPSILON)) {
            epsilon = (Parameter) xo.getChild(EPSILON).getChild(Parameter.class);
        } else {
            epsilon = new Parameter.Default("epsilon", 0.0, 0.0, Double.MAX_VALUE);
        }

        // --- ploidy ---
        Parameter ploidyFactors;
        if (xo.hasChildNamed(PLOIDY)) {
            ploidyFactors = (Parameter) xo.getChild(PLOIDY).getChild(Parameter.class);
        } else {
            ploidyFactors = new Parameter.Default(PLOIDY, intervalsList.size());
            for (int i = 0; i < intervalsList.size(); i++) {
                ploidyFactors.setParameterValue(i, 1.0);
            }
        }

        return new SquaredSplineCoalescentLikelihood(
                intervalsList, splines, epsilon, rateParameter, ploidyFactors);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new OrRule(
                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE),
                    new ElementRule(BigFastTreeIntervals.class, 1, Integer.MAX_VALUE)
            ),
            new ElementRule(SPLINES, new XMLSyntaxRule[]{
                    new ElementRule(SquaredCachedSplines.class)
            }),
            new ElementRule(RATE_PARAM, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),                                          // required
            new ElementRule(EPSILON, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),                                    // optional
            new ElementRule(PLOIDY, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),                                    // optional
    };

    @Override
    public String getParserDescription() {
        return "Coalescent likelihood for a squared-spline rate 1/N(t) = ε + (intercept + θ'b(t))². " +
               "The rateParameter=[coefficients,intercept] must be supplied explicitly so the same " +
               "object is shared by the HMC operator and any prior gradient wrapper.";
    }

    @Override
    public Class getReturnType() {
        return SquaredSplineCoalescentLikelihood.class;
    }
}
