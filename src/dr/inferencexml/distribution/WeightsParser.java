/*
 * WeightsParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.GriddedWeights;
import dr.inference.distribution.RandomField;
import dr.inference.distribution.Weights;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;

import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;

public class WeightsParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "weightProvider";
    private static final String RESCALE_BY_ROOT_HEIGHT = "rescaleByRootHeight";
    private static final String GRID_POINTS = "gridPoints";

    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        boolean rescaleByRootHeight = xo.getAttribute(RESCALE_BY_ROOT_HEIGHT, false);

        if (xo.hasChildNamed(GRID_POINTS)) {
            Parameter gridPoints = (Parameter) xo.getElementFirstChild(GRID_POINTS);
            return new GriddedWeights(tree, gridPoints, rescaleByRootHeight);
        } else {
            return new Weights(tree, rescaleByRootHeight);

        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(RESCALE_BY_ROOT_HEIGHT, true),
            new ElementRule(TreeModel.class),
            new ElementRule(GRID_POINTS, Parameter.class, "provide grid points", true)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Weights.class;
    }



}
