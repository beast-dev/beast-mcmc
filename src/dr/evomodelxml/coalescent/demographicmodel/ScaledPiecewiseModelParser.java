/*
 * ScaledPiecewiseModelParser.java
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

package dr.evomodelxml.coalescent.demographicmodel;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodel.PiecewisePopulationModel;
import dr.evomodel.coalescent.demographicmodel.ScaledPiecewiseModel;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a PiecewisePopulation.
 */
public class ScaledPiecewiseModelParser extends AbstractXMLObjectParser {

    public static final String PIECEWISE_POPULATION = "scaledPiecewisePopulation";
    public static final String EPOCH_SIZES = "populationSizes";
    public static final String TREE_MODEL = "populationTree";

    public String getParserName() {
        return PIECEWISE_POPULATION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(EPOCH_SIZES);
        Parameter epochSizes = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TREE_MODEL);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        boolean isLinear = xo.getBooleanAttribute("linear");

        return new ScaledPiecewiseModel(epochSizes, treeModel, isLinear, units);
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a piecewise population model";
    }

    public Class getReturnType() {
        return PiecewisePopulationModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(EPOCH_SIZES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(TREE_MODEL, new XMLSyntaxRule[]{new ElementRule(TreeModel.class)}),
            XMLUnits.SYNTAX_RULES[0],
            AttributeRule.newBooleanRule("linear")
    };

}
