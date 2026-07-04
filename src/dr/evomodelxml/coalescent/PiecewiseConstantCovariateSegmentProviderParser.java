/*
 * PiecewiseConstantCovariateSegmentProviderParser.java
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

import dr.evomodel.coalescent.timeline.PiecewiseConstantCovariateSegmentProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * XML parser for {@link PiecewiseConstantCovariateSegmentProvider}.
 *
 * Each segment k covers [breaks[k], breaks[k+1]).
 * The covariates parameter is a flat row-major vector: values for segment 0, then segment 1, ...
 * The {@code covariateDimension} attribute (default 1) says how many values per segment.
 *
 * Example (1 covariate, 3 segments over [0, 60000] ya):
 *
 * <pre>{@code
 * <piecewiseConstantCovariateSegmentProvider id="tempProxy">
 *     <breaks>
 *         <parameter value="0.0 20000.0 40000.0 60000.0"/>
 *     </breaks>
 *     <covariates covariateDimension="1">
 *         <parameter value="-1.5 -0.3 0.8"/>
 *     </covariates>
 * </piecewiseConstantCovariateSegmentProvider>
 * }</pre>
 *
 * @author Filippo Monti
 */
public class PiecewiseConstantCovariateSegmentProviderParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME        = "piecewiseConstantCovariateSegmentProvider";
    public static final String BREAKS             = "breaks";
    public static final String COVARIATES         = "covariates";
    public static final String COVARIATE_DIM_ATTR = "covariateDimension";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // --- breaks ---
        XMLObject breaksXO = xo.getChild(BREAKS);
        Parameter breaksParam = (Parameter) breaksXO.getChild(Parameter.class);
        int nBreaks = breaksParam.getDimension();
        if (nBreaks < 2) {
            throw new XMLParseException(PARSER_NAME + ": breaks must have at least 2 entries; got " + nBreaks);
        }
        double[] breaks = new double[nBreaks];
        for (int i = 0; i < nBreaks; i++) breaks[i] = breaksParam.getParameterValue(i);

        int nSegments = nBreaks - 1;

        // --- covariates ---
        XMLObject covXO = xo.getChild(COVARIATES);
        int covariateDim = covXO.getAttribute(COVARIATE_DIM_ATTR, 1);
        if (covariateDim < 1) {
            throw new XMLParseException(PARSER_NAME + ": covariateDimension must be >= 1; got " + covariateDim);
        }

        Parameter flatParam = (Parameter) covXO.getChild(Parameter.class);
        int expectedLen = nSegments * covariateDim;
        if (flatParam.getDimension() != expectedLen) {
            throw new XMLParseException(PARSER_NAME + ": covariates parameter has dimension " +
                    flatParam.getDimension() + " but expected nSegments * covariateDimension = " +
                    nSegments + " * " + covariateDim + " = " + expectedLen);
        }

        // reshape row-major: covariates[seg][dim]
        double[][] covariates = new double[nSegments][covariateDim];
        for (int seg = 0; seg < nSegments; seg++) {
            for (int d = 0; d < covariateDim; d++) {
                covariates[seg][d] = flatParam.getParameterValue(seg * covariateDim + d);
            }
        }

        return new PiecewiseConstantCovariateSegmentProvider(breaks, covariates);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BREAKS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(COVARIATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };

    @Override
    public String getParserDescription() {
        return "Piecewise-constant covariate segment provider. " +
               "Each segment [breaks[k], breaks[k+1]) carries a constant covariate vector. " +
               "The covariates parameter is flat row-major: covariateDimension values per segment.";
    }

    @Override
    public Class getReturnType() { return PiecewiseConstantCovariateSegmentProvider.class; }
}
