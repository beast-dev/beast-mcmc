/*
 * HierarchicalTransmissionDemographicModel.java
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

package dr.evomodel.transmission;

import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A demographic model for within patient evolution in a transmission history with a hierarchical parameterization.
 *
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id: TransmissionDemographicModel.java,v 1.7 2005/05/24 20:25:58 rambaut Exp $
 */
public class HierarchicalTransmissionDemographicModel extends TransmissionDemographicModel {

    //
    // Public stuff
    //

    public static String HIERARCHICAL_TRANSMISSION_MODEL = "hierarchicalTransmissionModel";

    public static String CONSTANT = TransmissionDemographicModel.CONSTANT;
    public static String EXPONENTIAL = TransmissionDemographicModel.EXPONENTIAL;
    public static String LOGISTIC = TransmissionDemographicModel.LOGISTIC;

    public static String POPULATION_SIZE = TransmissionDemographicModel.POPULATION_SIZE;
    public static String ANCESTRAL_PROPORTION = TransmissionDemographicModel.ANCESTRAL_PROPORTION;
    public static String GROWTH_RATE = TransmissionDemographicModel.GROWTH_RATE;
    public static String DOUBLING_TIME = TransmissionDemographicModel.DOUBLING_TIME;

    public static final String HOST_COUNT = "hostCount";

    public HierarchicalTransmissionDemographicModel(int model,
                                                    Parameter N0Parameter, Parameter N1Parameter,
                                                    Parameter growthRateParameter, Parameter doublingTimeParameter, Type units) {

        this(HIERARCHICAL_TRANSMISSION_MODEL, model, N0Parameter, N1Parameter, growthRateParameter, doublingTimeParameter, units);
    }

    public HierarchicalTransmissionDemographicModel(String name, int model,
                                                    Parameter N0Parameter, Parameter N1Parameter,
                                                    Parameter growthRateParameter, Parameter doublingTimeParameter, Type units) {

        super(name, model, N0Parameter, N1Parameter, growthRateParameter, doublingTimeParameter, units);
    }

    public void resizeHierarchicalParameters(int hostCount) {

        if (N0Parameter != null) {
            N0Parameter.setDimension(hostCount);
//            N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, hostCount));
        }

        if (N1Parameter != null) {
            N1Parameter.setDimension(hostCount);
//            N1Parameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        }

        if (growthRateParameter != null) {
            growthRateParameter.setDimension(hostCount);
//            growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, hostCount));
        }

        if (doublingTimeParameter != null) {
            doublingTimeParameter.setDimension(hostCount);
//            doublingTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, hostCount));
        }
    }

    protected int getIndexFromHost(int host) {
        return host - 1;
    }

    /**
     * Parses an element from an DOM document into a ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return HIERARCHICAL_TRANSMISSION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLUnits.Utils.getUnitsAttr(xo);

            int model = 0;

            Parameter N0Param = null;
            Parameter N1Param = null;
            Parameter rParam = null;
            Parameter dParam = null;

            if (xo.hasChildNamed(CONSTANT)) {

                XMLObject cxo = xo.getChild(CONSTANT);

                N0Param = (Parameter) cxo.getElementFirstChild(POPULATION_SIZE);
                model = 0;

            } else if (xo.hasChildNamed(EXPONENTIAL)) {

                XMLObject cxo = xo.getChild(EXPONENTIAL);

                N1Param = (Parameter) cxo.getElementFirstChild(ANCESTRAL_PROPORTION);
                if (cxo.hasChildNamed(GROWTH_RATE)) {
                    rParam = (Parameter) cxo.getElementFirstChild(GROWTH_RATE);
                } else {
                    dParam = (Parameter) cxo.getElementFirstChild(DOUBLING_TIME);
                }
                model = 1;

            } else if (xo.hasChildNamed(LOGISTIC)) {

                XMLObject cxo = xo.getChild(LOGISTIC);

                N0Param = (Parameter) cxo.getElementFirstChild(POPULATION_SIZE);
                N1Param = (Parameter) cxo.getElementFirstChild(ANCESTRAL_PROPORTION);

                if (cxo.hasChildNamed(GROWTH_RATE)) {
                    rParam = (Parameter) cxo.getElementFirstChild(GROWTH_RATE);
                } else {
                    dParam = (Parameter) cxo.getElementFirstChild(DOUBLING_TIME);
                }
                model = 2;
            }

            HierarchicalTransmissionDemographicModel demoModel =
                    new HierarchicalTransmissionDemographicModel(model, N0Param, N1Param, rParam, dParam, units);
            demoModel.resizeHierarchicalParameters(xo.getAttribute(HOST_COUNT, 1));

            return demoModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A SiteModel that has a gamma distributed rates across sites";
        }

        public Class getReturnType() {
            return HierarchicalTransmissionDemographicModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                XMLUnits.UNITS_RULE,
                AttributeRule.newIntegerRule(HOST_COUNT, false),
                new XORRule(
                        new ElementRule(CONSTANT,
                                new XMLSyntaxRule[]{
                                        new ElementRule(POPULATION_SIZE,
                                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                "This parameter represents the carrying capacity (maximum population size). " +
                                                        "If the shape is very large then the current day population size will be very close to the carrying capacity."),
                                }
                        ),
                        new XORRule(
                                new ElementRule(EXPONENTIAL,
                                        new XMLSyntaxRule[]{
                                                new XORRule(
                                                        new ElementRule(GROWTH_RATE,
                                                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                                "This parameter determines the rate of growth during the exponential phase. See exponentialGrowth for details."),
                                                        new ElementRule(DOUBLING_TIME,
                                                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                                "This parameter determines the doubling time at peak growth rate.")),
                                                new ElementRule(ANCESTRAL_PROPORTION,
                                                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                        "This parameter determines the populaation size at transmission.")
                                        }
                                ),
                                new ElementRule(LOGISTIC,
                                        new XMLSyntaxRule[]{
                                                new ElementRule(POPULATION_SIZE,
                                                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                        "This parameter represents the carrying capacity (maximum population size). " +
                                                                "If the shape is very large then the current day population size will be very close to the carrying capacity."),
                                                new XORRule(
                                                        new ElementRule(GROWTH_RATE,
                                                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                                "This parameter determines the rate of growth during the exponential phase. See exponentialGrowth for details."),
                                                        new ElementRule(DOUBLING_TIME,
                                                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                                "This parameter determines the doubling time at peak growth rate.")),
                                                new ElementRule(ANCESTRAL_PROPORTION,
                                                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                                        "This parameter determines the populaation size at transmission.")
                                        }
                                )
                        )
                )
        };

    };

}
