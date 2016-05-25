/*
 * TransmissionDemographicModel.java
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

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evomodel.coalescent.DemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A demographic model for within patient evolution in a transmission history.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TransmissionDemographicModel.java,v 1.7 2005/05/24 20:25:58 rambaut Exp $
 */
public class TransmissionDemographicModel extends DemographicModel {

    //
    // Public stuff
    //

    public static String TRANSMISSION_MODEL = "transmissionModel";

    public static String CONSTANT = "constant";
    public static String EXPONENTIAL = "exponential";
    public static String LOGISTIC = "logistic";

    public static String POPULATION_SIZE = "populationSize";
    public static String ANCESTRAL_PROPORTION = "ancestralProportion";
    public static String GROWTH_RATE = "growthRate";
    public static String DOUBLING_TIME = "doublingTime";

    /**
     * Construct demographic model with default settings
     */
    public TransmissionDemographicModel(int model,
                                        Parameter N0Parameter, Parameter N1Parameter,
                                        Parameter growthRateParameter, Parameter doublingTimeParameter, Type units) {

        this(TRANSMISSION_MODEL, model, N0Parameter, N1Parameter, growthRateParameter, doublingTimeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public TransmissionDemographicModel(String name, int model,
                                        Parameter N0Parameter, Parameter N1Parameter,
                                        Parameter growthRateParameter, Parameter doublingTimeParameter, Type units) {

        super(name);

        this.model = model;
        if (model == 0) {
            hostDemographic = new ConstantPopulation(units);
        } else if (model == 1) {
            hostDemographic = new TransmissionExponentialGrowth(units);
        } else if (model == 2) {
            hostDemographic = new TransmissionLogisticGrowth(units);
        }

        if (N0Parameter != null) {
            this.N0Parameter = N0Parameter;
            addVariable(N0Parameter);
            N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        if (N1Parameter != null) {
            this.N1Parameter = N1Parameter;
            addVariable(N1Parameter);
            N1Parameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        }

        if (growthRateParameter != null) {
            this.growthRateParameter = growthRateParameter;
            addVariable(growthRateParameter);
            growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        if (doublingTimeParameter != null) {
            this.doublingTimeParameter = doublingTimeParameter;
            addVariable(doublingTimeParameter);
            doublingTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        setUnits(units);
    }

    public DemographicFunction getDemographicFunction() {
        throw new RuntimeException("getDemographicFunction not used for TransmissionDemographicModel");
    }

    protected int getIndexFromHost(int host) {
        return 0;
    }

    // Attempting to get the transmission model to be host-specific
    public DemographicFunction getDemographicFunction(double transmissionTime, double donorSize, int host) {

        final int index = getIndexFromHost(host);

        if (model == 0) {
            // constant
            double N0 = N0Parameter.getParameterValue(index);
            ((ConstantPopulation) hostDemographic).setN0(N0);

        } else if (model == 1) {
            // exponential
            double N1 = N1Parameter.getParameterValue(index);
            ((TransmissionDemographicFunction) hostDemographic).setTransmissionTime(transmissionTime);
            ((TransmissionDemographicFunction) hostDemographic).setDonorSize(donorSize);
            ((TransmissionDemographicFunction) hostDemographic).setBottleNeckProportion(N1);
            if (growthRateParameter != null) {
                double r = growthRateParameter.getParameterValue(index);
                ((TransmissionExponentialGrowth) hostDemographic).setGrowthRate(r);
            } else {
                double d = doublingTimeParameter.getParameterValue(index);
                ((TransmissionExponentialGrowth) hostDemographic).setDoublingTime(d);
            }

        } else if (model == 2) {
            // logistic
            ((TransmissionDemographicFunction) hostDemographic).setTransmissionTime(transmissionTime);
            ((TransmissionDemographicFunction) hostDemographic).setDonorSize(donorSize);

            double N0 = N0Parameter.getParameterValue(index);
            ((TransmissionLogisticGrowth) hostDemographic).setN0(N0);
            double N1 = N1Parameter.getParameterValue(index);
            ((TransmissionDemographicFunction) hostDemographic).setBottleNeckProportion(N1);
            if (growthRateParameter != null) {
                double r = growthRateParameter.getParameterValue(index);
                ((TransmissionLogisticGrowth) hostDemographic).setGrowthRate(r);
            } else {
                double d = doublingTimeParameter.getParameterValue(index);
                ((TransmissionLogisticGrowth) hostDemographic).setDoublingTime(d);
            }

        }

        return hostDemographic;
    }

    /**
     * Parses an element from an DOM document into a ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLUnits.Utils.getUnitsAttr(xo);

            int model = 0;

            Parameter N0Param = null;
            Parameter N1Param = null;
            Parameter rParam = null;
            Parameter dParam = null;

            if (xo.hasChildNamed(CONSTANT)) {

                XMLObject cxo = (XMLObject) xo.getChild(CONSTANT);

                N0Param = (Parameter) cxo.getElementFirstChild(POPULATION_SIZE);
                model = 0;

            } else if (xo.hasChildNamed(EXPONENTIAL)) {

                XMLObject cxo = (XMLObject) xo.getChild(EXPONENTIAL);

                N1Param = (Parameter) cxo.getElementFirstChild(ANCESTRAL_PROPORTION);
                if (cxo.hasChildNamed(GROWTH_RATE)) {
                    rParam = (Parameter) cxo.getElementFirstChild(GROWTH_RATE);
                } else {
                    dParam = (Parameter) cxo.getElementFirstChild(DOUBLING_TIME);
                }
                model = 1;

            } else if (xo.hasChildNamed(LOGISTIC)) {

                XMLObject cxo = (XMLObject) xo.getChild(LOGISTIC);

                N0Param = (Parameter) cxo.getElementFirstChild(POPULATION_SIZE);
                N1Param = (Parameter) cxo.getElementFirstChild(ANCESTRAL_PROPORTION);

                if (cxo.hasChildNamed(GROWTH_RATE)) {
                    rParam = (Parameter) cxo.getElementFirstChild(GROWTH_RATE);
                } else {
                    dParam = (Parameter) cxo.getElementFirstChild(DOUBLING_TIME);
                }
                model = 2;
            }

            return new TransmissionDemographicModel(model, N0Param, N1Param, rParam, dParam, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A SiteModel that has a gamma distributed rates across sites";
        }

        public Class getReturnType() {
            return TransmissionDemographicModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                XMLUnits.UNITS_RULE,
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

    private abstract class TransmissionDemographicFunction extends DemographicFunction.Abstract {
        public TransmissionDemographicFunction(Type units) {
            super(units);
        }

        public void setDonorSize(double donorSize) {
            this.donorSize = donorSize;
        }

        public void setTransmissionTime(double transmissionTime) {
            this.transmissionTime = transmissionTime;
        }

        public abstract void setBottleNeckProportion(double prop);

        protected double transmissionTime;
        protected double donorSize;
    }

    private class TransmissionExponentialGrowth extends TransmissionDemographicFunction {

        public TransmissionExponentialGrowth(Type units) {
            super(units);
        }

        public void setGrowthRate(double r) {
            this.r = r;
        }

        public void setDoublingTime(double doublingTime) {
            setGrowthRate(Math.log(2) / doublingTime);
        }

        public void setBottleNeckProportion(double prop) {
            this.N1 = prop * donorSize;
        }

        // Implementation of abstract methods

        public double getDemographic(double t) {

            return N1 * Math.exp(-r * (t - transmissionTime));
        }

        public double getIntegral(double start, double finish) {

            double t1 = start - transmissionTime;
            double t2 = finish - transmissionTime;

            double integral = ((Math.exp(t2 * r) - 1.0) / N1 / r) - ((Math.exp(t1 * r) - 1.0) / N1 / r);

            return integral;
        }

        public double getIntensity(double t) {
            throw new RuntimeException("Function not used");
        }

        public double getInverseIntensity(double x) {
            throw new RuntimeException("Function not used");
        }

        public int getNumArguments() {
            throw new RuntimeException("Function not used");
        }

        public String getArgumentName(int n) {
            throw new RuntimeException("Function not used");
        }

        public double getArgument(int n) {
            throw new RuntimeException("Function not used");
        }

        public void setArgument(int n, double value) {
            throw new RuntimeException("Function not used");
        }

        public double getLowerBound(int n) {
            throw new RuntimeException("Function not used");
        }

        public double getUpperBound(int n) {
            throw new RuntimeException("Function not used");
        }

        public DemographicFunction getCopy() {
            TransmissionExponentialGrowth t = new TransmissionExponentialGrowth(getUnits());
            t.r = r;
            t.N1 = N1;
            t.transmissionTime = transmissionTime;
            t.donorSize = donorSize;

            return t;
        }

        private double r;
        private double N1;
    }

    private class TransmissionLogisticGrowth extends TransmissionDemographicFunction {

        public TransmissionLogisticGrowth(Type units) {
            super(units);
        }

        public void setGrowthRate(double r) {
            this.r = r;
        }

        public void setDoublingTime(double doublingTime) {
            setGrowthRate(Math.log(2) / doublingTime);
        }

        public void setBottleNeckProportion(double prop) {
            this.N1 = prop * donorSize;
        }

        public void setN0(double N0) {
            this.N0 = N0;
        }

        public double getN0() {
            return N0;
        }

        // Implementation of abstract methods

        public double getDemographic(double t) {

            double common = Math.exp(-r * (t - transmissionTime - (Math.log(1.0 / ((N0 / N1) - 1.0)) / r)));
            return (N0 * common) / (1.0 + common);
        }

        public double getIntegral(double start, double finish) {

            double g = finish - start;
            double t = start - transmissionTime - (Math.log(1.0 / ((N0 / N1) - 1.0)) / r);

            double eMinusRT = Math.exp(-r * t);
            double eMinusRG = Math.exp(-r * g);

            double integral = (g / N0) + ((1 - eMinusRG) / (N0 * eMinusRT * r * eMinusRG));

//			double integral2 = getNumericalIntegral(start, finish);
//			if (Math.abs(integral - integral2) > 1E-8) {
//				System.err.println("Numerical integration failed: " + integral + " (truth=" + integral2 + ")");
//				//throw new RuntimeException("Numerical integration failed: " + integral + " (truth=" + integral2 + ")");
//			}

            return integral;
        }

        public double getIntensity(double t) {
            throw new RuntimeException("Function not used");
        }

        public double getInverseIntensity(double x) {
            throw new RuntimeException("Function not used");
        }

        public int getNumArguments() {
            throw new RuntimeException("Function not used");
        }

        public String getArgumentName(int n) {
            throw new RuntimeException("Function not used");
        }

        public double getArgument(int n) {
            throw new RuntimeException("Function not used");
        }

        public void setArgument(int n, double value) {
            throw new RuntimeException("Function not used");
        }

        public double getLowerBound(int n) {
            throw new RuntimeException("Function not used");
        }

        public double getUpperBound(int n) {
            throw new RuntimeException("Function not used");
        }

        public DemographicFunction getCopy() {
            TransmissionLogisticGrowth t = new TransmissionLogisticGrowth(getUnits());
            t.N0 = N0;
            t.r = r;
            t.N1 = N1;
            t.transmissionTime = transmissionTime;
            t.donorSize = donorSize;

            return t;
        }

        private double N0;
        private double r;
        private double N1;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter N1Parameter = null;
    Parameter growthRateParameter = null;
    Parameter doublingTimeParameter = null;
    DemographicFunction hostDemographic = null;

    int model = 0;
}
