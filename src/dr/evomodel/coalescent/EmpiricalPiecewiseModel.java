/*
 * EmpiricalPiecewiseModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.EmpiricalPiecewiseConstant;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: EmpiricalPiecewiseModel.java,v 1.4 2005/04/11 11:24:39 alexei Exp $
 */
public class EmpiricalPiecewiseModel extends DemographicModel {
    //
    // Public stuff
    //

    public static final String EMPIRICAL_PIECEWISE = "empiricalPiecewise";
    public static final String INTERVAL_WIDTHS = "intervalWidths";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String TAU = "generationLength";
    public static final String THRESHOLD = "threshold";
    public static final String LAG = "lag";

    /**
     * Construct demographic model with default settings
     */
    public EmpiricalPiecewiseModel(double[] intervalWidths, Parameter populationSizesParameter, Parameter tauParameter, Parameter bParameter, Parameter lagParameter, Type units) {

        this(EMPIRICAL_PIECEWISE, intervalWidths, populationSizesParameter, tauParameter, bParameter, lagParameter, units);

    }

    /**
     * Construct demographic model with default settings
     */
    public EmpiricalPiecewiseModel(String name, double[] intervalWidths, Parameter populationSizesParameter, Parameter tauParameter, Parameter bParameter, Parameter lagParameter, Type units) {

        super(name);

        //System.out.println("intervalWidths.length=" + intervalWidths.length);
        //System.out.println("populationSizes.dimension=" + populationSizesParameter.getDimension());

        if (intervalWidths.length == 1) {
            double[] newIntervalWidths = new double[populationSizesParameter.getDimension() - 1];
            for (int i = 0; i < newIntervalWidths.length; i++) {
                newIntervalWidths[i] = intervalWidths[0];
            }
            intervalWidths = newIntervalWidths;
        }
        //System.out.println("new intervalWidths.length=" + intervalWidths.length);

        if (populationSizesParameter.getDimension() != (intervalWidths.length + 1)) {
            throw new IllegalArgumentException(
                    "interval widths array must have either 1 or " + (populationSizesParameter.getDimension() - 1) +
                            " elements, but instead it has " + intervalWidths.length + "."
            );
        }

        this.tauParameter = tauParameter;
        this.lagParameter = lagParameter;
        this.bParameter = bParameter;

        this.intervalWidths = intervalWidths;
        this.populationSizesParameter = populationSizesParameter;

        addVariable(tauParameter);
        addVariable(lagParameter);
        addVariable(bParameter);
        addVariable(populationSizesParameter);
        tauParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, tauParameter.getDimension()));
        lagParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, lagParameter.getDimension()));
        bParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, bParameter.getDimension()));
        populationSizesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, populationSizesParameter.getDimension()));

        setUnits(units);

        piecewiseFunction = new EmpiricalPiecewiseConstant(intervalWidths, calculatePopSizes(),
                lagParameter.getParameterValue(0), units);
    }

    /**
     *
     */
    public DemographicFunction getDemographicFunction() {
        piecewiseFunction.setLag(lagParameter.getParameterValue(0));
        piecewiseFunction.setPopulationSizes(calculatePopSizes());

        return piecewiseFunction;
    }

    private double[] calculatePopSizes() {
        double m = tauParameter.getParameterValue(0);
        double c = bParameter.getParameterValue(0);

        double[] popSizes = new double[populationSizesParameter.getDimension()];
        for (int i = 0; i < popSizes.length; i++) {
            popSizes[i] = m * populationSizesParameter.getParameterValue(i) + c;
        }
        return popSizes;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

//	protected void handleVariableChangedEvent(Parameter parameter, int index) {
//
//		// no intermediates need to be recalculated...
//	}
//

    // todo: why override?

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    /**
     * Parses an element from an DOM document into a PiecewisePopulation.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EMPIRICAL_PIECEWISE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLUnits.Utils.getUnitsAttr(xo);


            XMLObject cxo = xo.getChild(INTERVAL_WIDTHS);
            double[] intervalWidths = cxo.getDoubleArrayAttribute("values");

            cxo = xo.getChild(POPULATION_SIZES);
            Parameter popSizes = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(TAU);
            Parameter scaleParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(THRESHOLD);
            Parameter bParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(LAG);
            Parameter offsetParam = (Parameter) cxo.getChild(Parameter.class);

            return new EmpiricalPiecewiseModel(intervalWidths, popSizes, scaleParam, bParam, offsetParam, units);
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
                new ElementRule(INTERVAL_WIDTHS,
                        new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule("values", false),}),
                XMLUnits.SYNTAX_RULES[0],
                new ElementRule(POPULATION_SIZES,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                        "The effective population sizes of each interval."),
                new ElementRule(TAU,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                        "The scale factor."),
                new ElementRule(THRESHOLD,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                        "The threshold before counts occur."),
                new ElementRule(LAG,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                        "The lag between actual population sizes and genetic diversity.")
        };
    };


    //
    // protected stuff
    //

    Parameter tauParameter;
    Parameter lagParameter;
    Parameter bParameter;
    Parameter populationSizesParameter;
    double[] intervalWidths;
    EmpiricalPiecewiseConstant piecewiseFunction = null;
}
