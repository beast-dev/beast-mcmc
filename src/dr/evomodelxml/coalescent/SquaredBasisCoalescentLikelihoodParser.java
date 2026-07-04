/*
 * SquaredBasisCoalescentLikelihoodParser.java
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
import dr.evomodel.coalescent.SquaredBasisCoalescentLikelihood;
import dr.evomodel.coalescent.basis.BSplineBasisExpansion;
import dr.evomodel.coalescent.basis.CovariateAugmentedBasisExpansion;
import dr.evomodel.coalescent.basis.CovariateMode;
import dr.evomodel.coalescent.timeline.CoalescentSegmentProvider;
import dr.evomodel.coalescent.timeline.NoCovariateSegmentProvider;
import dr.evomodel.coalescent.timeline.PiecewiseConstantCovariateSegmentProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import dr.xml.*;

import java.util.Locale;

import java.util.ArrayList;
import java.util.List;

/**
 * XML parser for SquaredBasisCoalescentLikelihood.
 *
 * <h3>No-covariate (splines only)</h3>
 *
 * The {@code <rateParameter>} must be a CompoundParameter in CANONICAL layout
 * {@code [intercept, coefficients...]} (intercept first).
 * The parameter is shared with the HMC operator and any prior gradient wrapper.
 *
 * <pre>{@code
 * <compoundParameter id="squaredBasisRate">
 *     <parameter idref="spline.intercept"/>
 *     <parameter idref="spline.coefficients"/>
 * </compoundParameter>
 *
 * <squaredBasisCoalescentLikelihood id="coalescent">
 *     <treeModel idref="treeModel"/>
 *     <splines><squaredCachedSplines idref="splineRates"/></splines>
 *     <rateParameter><compoundParameter idref="squaredBasisRate"/></rateParameter>
 *     <epsilon><parameter idref="epsilon"/></epsilon>
 * </squaredBasisCoalescentLikelihood>
 * }</pre>
 *
 * <h3>With covariates (ADDITIVE)</h3>
 *
 * Add an optional {@code <covariateModel>} block. The {@code <rateParameter>} still contains
 * ONLY the spline part ({@code [intercept, coefficients...]}). Covariate coefficients are
 * declared separately. The likelihood's parameter (accessible to HMC) is the INTERNAL compound
 * {@code [rateParameter, covariateCoefficients]}, created automatically.
 *
 * For HMC, define a matching compound in XML that shares the same sub-parameter objects:
 *
 * <pre>{@code
 * <compoundParameter id="squaredBasisRate">
 *     <parameter idref="spline.intercept"/>
 *     <parameter idref="spline.coefficients"/>
 * </compoundParameter>
 *
 * <parameter id="cov.coeff" value="0.0"/>
 *
 * <!-- fullBasisRate shares sub-params with the likelihood's internal compound -->
 * <compoundParameter id="fullBasisRate">
 *     <compoundParameter idref="squaredBasisRate"/>
 *     <parameter idref="cov.coeff"/>
 * </compoundParameter>
 *
 * <squaredBasisCoalescentLikelihood id="coalescent">
 *     <treeModel idref="treeModel"/>
 *     <splines><squaredCachedSplines idref="splineRates"/></splines>
 *     <rateParameter><compoundParameter idref="squaredBasisRate"/></rateParameter>
 *     <epsilon><parameter idref="epsilon"/></epsilon>
 *     <covariateModel>
 *         <piecewiseConstantCovariateSegmentProvider idref="tempProxy"/>
 *         <covariateCoefficients><parameter idref="cov.coeff"/></covariateCoefficients>
 *     </covariateModel>
 * </squaredBasisCoalescentLikelihood>
 *
 * <!-- HMC uses fullBasisRate; events propagate through shared sub-parameters -->
 * <hamiltonianMonteCarloOperator ...>
 *     <squaredBasisCoalescentLikelihood idref="coalescent"/>
 *     <compoundParameter idref="fullBasisRate"/>
 * </hamiltonianMonteCarloOperator>
 * }</pre>
 *
 * @author Filippo Monti
 */
public class SquaredBasisCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME      = "squaredBasisCoalescentLikelihood";
    public static final String SPLINES          = "splines";
    public static final String RATE_PARAM       = "rateParameter";
    public static final String EPSILON          = "epsilon";
    public static final String PLOIDY           = "ploidy";
    public static final String COVARIATE_MODEL  = "covariateModel";
    public static final String COV_COEFFICIENTS = "covariateCoefficients";
    public static final String COVARIATE_MODE   = "covariateMode";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // --- splines ---
        XMLObject splinesXO = xo.getChild(SPLINES);
        SquaredCachedSplines splines =
                (SquaredCachedSplines) splinesXO.getChild(SquaredCachedSplines.class);

        // --- rateParameter: canonical spline compound [intercept, coefficients...] ---
        XMLObject rateParamXO = xo.getChild(RATE_PARAM);
        Parameter rateParameter = (Parameter) rateParamXO.getChild(Parameter.class);

        // Validate spline dimension
        int splineDim = splines.getCoefficientDim() + 1;  // +1 for intercept
        if (rateParameter.getDimension() != splineDim) {
            throw new XMLParseException(PARSER_NAME + ": rateParameter dimension " +
                    rateParameter.getDimension() + " must equal spline coefficients + 1 = " + splineDim);
        }

        // Validate canonical layout: compound must wrap splines.getCoefficients() and getIntercept()
        if (rateParameter instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) rateParameter;
            boolean coefFound = false, interceptFound = false;
            for (int k = 0; k < cp.getParameterCount(); k++) {
                Parameter sub = cp.getParameter(k);
                if (sub == splines.getCoefficients()) coefFound      = true;
                if (sub == splines.getIntercept())    interceptFound = true;
            }
            if (!coefFound || !interceptFound) {
                throw new XMLParseException(PARSER_NAME + ": rateParameter must wrap " +
                        "splines.getCoefficients() and splines.getIntercept()");
            }
        }

        // --- trees ---
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
            for (int i = 0; i < intervalsList.size(); i++) ploidyFactors.setParameterValue(i, 1.0);
        }

        // --- covariateModel (optional) ---
        BSplineBasisExpansion bSplineBasis = new BSplineBasisExpansion(splines, rateParameter, true);

        if (xo.hasChildNamed(COVARIATE_MODEL)) {
            XMLObject covModel = xo.getChild(COVARIATE_MODEL);

            PiecewiseConstantCovariateSegmentProvider segProvider =
                    (PiecewiseConstantCovariateSegmentProvider)
                            covModel.getChild(PiecewiseConstantCovariateSegmentProvider.class);

            XMLObject covCoeffXO = covModel.getChild(COV_COEFFICIENTS);
            Parameter covCoeff = (Parameter) covCoeffXO.getChild(Parameter.class);

            // Parse optional covariateMode attribute (default: ADDITIVE).
            String modeStr = covModel.getAttribute(COVARIATE_MODE, "additive");
            CovariateMode covMode;
            try {
                covMode = CovariateMode.valueOf(modeStr.toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(PARSER_NAME + ": unknown covariateMode '" + modeStr +
                        "'; expected 'additive' or 'multiplicative'");
            }

            int covariateDim = segProvider.getCovariateDimension();
            int expectedCoeffDim = (covMode == CovariateMode.ADDITIVE)
                    ? covariateDim
                    : covariateDim * splineDim;
            if (covCoeff.getDimension() != expectedCoeffDim) {
                throw new XMLParseException(PARSER_NAME + ": covariateCoefficients dimension " +
                        covCoeff.getDimension() + " must equal " + expectedCoeffDim +
                        " for mode=" + covMode);
            }

            CovariateAugmentedBasisExpansion augBasis = new CovariateAugmentedBasisExpansion(
                    bSplineBasis, covCoeff, covariateDim, covMode);

            return new SquaredBasisCoalescentLikelihood(
                    intervalsList, augBasis, segProvider, epsilon, ploidyFactors);
        }

        // No covariates
        return new SquaredBasisCoalescentLikelihood(
                intervalsList, bSplineBasis, NoCovariateSegmentProvider.INSTANCE,
                epsilon, ploidyFactors);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

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
            }),
            new ElementRule(EPSILON, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(PLOIDY, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(COVARIATE_MODEL, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(COVARIATE_MODE, true),
                    new ElementRule(PiecewiseConstantCovariateSegmentProvider.class),
                    new ElementRule(COV_COEFFICIENTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
            }, true),
    };

    @Override
    public String getParserDescription() {
        return "Coalescent likelihood for a squared-basis rate 1/N(t) = ε + (γ'ψ(t))². " +
               "rateParameter=[intercept, coefficients...] is in canonical layout (intercept first). " +
               "Optional covariateModel block adds piecewise-constant additive covariates via " +
               "CovariateAugmentedBasisExpansion (ADDITIVE mode).";
    }

    @Override
    public Class getReturnType() { return SquaredBasisCoalescentLikelihood.class; }
}
