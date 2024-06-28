/*
 * TimeVaryingBranchRateModelParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Parameter;
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel.FunctionalForm.Type;
import dr.xml.*;

import static dr.evomodelxml.coalescent.smooth.SmoothSkygridLikelihoodParser.getGridPoints;

public class TimeVaryingBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeVaryingRates";

    private static final String RATES = StrictClockBranchRatesParser.RATE;
    private static final String GRID_POINTS = GMRFSkyrideLikelihoodParser.GRID_POINTS;
    private static final String NUM_GRID_POINTS = GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS;
    private static final String CUT_OFF = GMRFSkyrideLikelihoodParser.CUT_OFF;

    private static final String FUNCTIONAL_FORM = "functionalForm";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        Parameter rates = (Parameter) xo.getElementFirstChild(RATES);
        Parameter gridPoints = getGridPoints(xo);

        if (rates.getDimension() != gridPoints.getDimension() + 1) {
            throw new XMLParseException("Rates dimension != gridPoints dimension + 1");
        }

        Type type = Type.parse(xo.getAttribute(FUNCTIONAL_FORM, Type.PIECEWISE_CONSTANT.getName()));

        return new TimeVaryingBranchRateModel(type, tree, rates, gridPoints);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TimeVaryingBranchRateModelParser.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new AndRule(
                            new ElementRule(CUT_OFF, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                            new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
                    )
            ),
            AttributeRule.newStringRule(FUNCTIONAL_FORM, true),
    };
}
