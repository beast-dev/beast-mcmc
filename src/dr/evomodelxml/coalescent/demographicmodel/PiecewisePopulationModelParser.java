/*
 * PiecewisePopulationModelParser.java
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

package dr.evomodelxml.coalescent.demographicmodel;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodel.PiecewisePopulationModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a PiecewisePopulation.
 */
public class PiecewisePopulationModelParser extends AbstractXMLObjectParser {

    public static final String PIECEWISE_POPULATION = "piecewisePopulation";
    public static final String EPOCH_SIZES = "epochSizes";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_RATES = "growthRates";
    public static final String EPOCH_WIDTHS = "epochWidths";
    public static final String GRID_POINTS = "gridPoints";

    public static final String WIDTHS = "widths";
    public static final String LINEAR = "linear";

    public String getParserName() {
        return PIECEWISE_POPULATION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        double[] epochWidths;

        if (xo.hasChildNamed(EPOCH_WIDTHS)) {
            XMLObject obj = xo.getChild(EPOCH_WIDTHS);
            epochWidths = obj.getDoubleArrayAttribute(WIDTHS);
        } else if (xo.hasChildNamed(GRID_POINTS)) {
            Parameter parameter = (Parameter) xo.getElementFirstChild(GRID_POINTS);
            double[] gridPoints = parameter.getParameterValues();
            epochWidths = computeWidths(gridPoints);
        } else {
            throw new XMLParseException("PiecewisePopulationModel must have an epochWidths or GridPoints element.");
        }

//        XMLObject obj = xo.getChild(EPOCH_WIDTHS);
//        double[] epochWidths = obj.getDoubleArrayAttribute(WIDTHS);

        if (xo.hasChildNamed(EPOCH_SIZES)) {
            Parameter epochSizes = (Parameter) xo.getElementFirstChild(EPOCH_SIZES);
            boolean isLinear = false;
            if (xo.hasAttribute(LINEAR)) {
                isLinear = xo.getBooleanAttribute(LINEAR);
            }
            return new PiecewisePopulationModel(PIECEWISE_POPULATION, epochSizes, epochWidths, isLinear, units);
        } else {
            Parameter populationSize = (Parameter) xo.getElementFirstChild(POPULATION_SIZE);
            if (xo.hasChildNamed(GROWTH_RATES)) {
                Parameter growthRates = (Parameter) xo.getElementFirstChild(GROWTH_RATES);
                return new PiecewisePopulationModel(PIECEWISE_POPULATION, populationSize, growthRates, epochWidths, units);
            } else {
                return new PiecewisePopulationModel(PIECEWISE_POPULATION, populationSize, epochWidths, false, units);

            }
        }
    }

    private double[] computeWidths(double[] gridPoints) throws XMLParseException {
        double epochWidths[] = new double[gridPoints.length];
        epochWidths[0] = gridPoints[0]; // TODO change this in case t_0 \neq 0
        for(int i = 1; i < gridPoints.length; i++) {
            if(gridPoints[i] <= gridPoints[i-1]) {
                throw new XMLParseException("Grid points must be in increasing order.");
            } else {
                epochWidths[i] = gridPoints[i] - gridPoints[i-1];
            }
        }
        return epochWidths;
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

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(EPOCH_SIZES,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new AndRule(
                            new ElementRule(POPULATION_SIZE,
//                                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                                    new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Parameter.class),
                                                    new ElementRule(Transform.ParsedTransform.class))}),
                            new ElementRule(GROWTH_RATES,
                                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
                    )
            ),
            new XORRule(
                    new ElementRule(EPOCH_WIDTHS,
                            new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule(WIDTHS)}),
                    new ElementRule(GRID_POINTS,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
            ),
            AttributeRule.newBooleanRule(LINEAR, true)
    };
}
