/*
 * MultivariateNormalIndependenceSampler.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inference.regression.SelfControlledCaseSeries;
import dr.math.MathUtils;
import dr.math.Poisson;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

import java.util.HashSet;
import java.util.Set;

/**
 * An independent normal distribution sampler to propose new (independent) values from a provided normal distribution model.
 *
 * @author Marc Suchard
 * @author Guy Baele
 */
public class MultivariateNormalIndependenceSampler extends AbstractCoercableOperator {

    public static final String OPERATOR_NAME = "multivariateNormalIndependenceSampler";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String SET_SIZE_MEAN = "setSizeMean";

    private double scaleFactor;
    private final Parameter parameter;
    private final int dim;
    private double setSizeMean;
    private final SelfControlledCaseSeries sccs;

    public MultivariateNormalIndependenceSampler(Parameter parameter,
                                                 SelfControlledCaseSeries sccs,
                                                 double setSizeMean,
                                                 double weight, double scaleFactor, CoercionMode mode) {
        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        setWeight(weight);
        dim = parameter.getDimension();
        setWeight(weight);
        this.sccs = sccs;
        this.setSizeMean = setSizeMean;
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return "independentNormalDistribution(" + parameter.getVariableName() + ")";
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public double doOperation() {

        double[] mean = sccs.getMode();
        double[] currentValue = parameter.getParameterValues();
        double[] newValue = new double[dim];

        Set<Integer> updateSet = new HashSet<Integer>();

        if (setSizeMean != -1.0) {
            final int listLength = Poisson.nextPoisson(setSizeMean);
            while (updateSet.size() <  listLength) {
                int newInt = MathUtils.nextInt(parameter.getDimension());
                if (!updateSet.contains(newInt)) {
                    updateSet.add(newInt);
                }
            }
        } else {
            for (int i = 0; i < dim; ++i) {
                updateSet.add(i);
            }
        }

        double logq = 0;
        for (Integer i : updateSet) {
            newValue[i] = mean[i] + scaleFactor * MathUtils.nextGaussian();
            if (UPDATE_ALL) {
                parameter.setParameterValueQuietly(i, newValue[i]);
            } else {
                parameter.setParameterValue(i, newValue[i]);
            }

            logq += (NormalDistribution.logPdf(currentValue[i], mean[i], scaleFactor) -
                    NormalDistribution.logPdf(newValue[i], mean[i], scaleFactor));
        }

//        for (Integer i : updateSet) {
//            parameter.setParameterValueQuietly(i, newValue[i]);
//        }

        if (UPDATE_ALL) {
            parameter.setParameterValueNotifyChangedAll(0, parameter.getParameterValue(0));
        }

        return logq;
    }

    private static final boolean UPDATE_ALL = false;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0) {
                throw new XMLParseException("scaleFactor must be greater than 0.0");
            }

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            SelfControlledCaseSeries sccs = (SelfControlledCaseSeries) xo.getChild(SelfControlledCaseSeries.class);

            double setSizeMean = xo.getAttribute(SET_SIZE_MEAN, -1.0);

            return new MultivariateNormalIndependenceSampler(parameter, sccs, setSizeMean, weight, scaleFactor, mode);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newDoubleRule(SET_SIZE_MEAN, true),
                new ElementRule(SelfControlledCaseSeries.class),
                new ElementRule(Parameter.class),
        };

        public String getParserDescription() {
            return "This element returns an independence sampler from a provided normal distribution model.";
        }

        public Class getReturnType() {
            return MultivariateNormalIndependenceSampler.class;
        }

    };

    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

}
