/*
 * ModeFindOperator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.multidimensionalscaling.mm.MMAlgorithm;
import dr.inference.multidimensionalscaling.mm.MultiDimensionalScalingMM;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

import javax.sql.rowset.serial.SerialRef;


/**
 * @author Marc Suchard
 */
public class ModeFindOperator extends AbstractCoercableOperator {

    public static final String OPERATOR = "modeFindOperator";
    public static final String MAX_TIMES = "maxTimes";

    private double scaleFactor;

    private final int maxTimes;
    private final int maxSteps;
    private int executeTimes = 0;
    private final MultiDimensionalScalingMM mm;

    public ModeFindOperator(MultiDimensionalScalingMM mm, int maxTimes, double weight, CoercionMode mode) {
        this(mm, maxTimes, weight, 1000, mode);
    }

    public ModeFindOperator(MultiDimensionalScalingMM mm, int maxTimes, double weight,
                            int maxModeSteps, CoercionMode mode) {
        super(mode);
        setWeight(weight);

        this.maxTimes = maxTimes;
        this.mm = mm;
        this.maxSteps = maxModeSteps;
    }

    public double doOperation() throws OperatorFailedException {

        if (executeTimes < maxTimes) {

            MatrixParameterInterface parameter = mm.getLikelihood().getMatrixParameter();

//            double[] original = null;

//            if (mode == CoercionMode.COERCION_ON) {
//                original = parameter.getParameterValues();
//            }

            System.err.println("START");
            mm.run(maxSteps);
            System.err.println("END");
            ++executeTimes;

            double logHR = Double.POSITIVE_INFINITY; // Always accept, breaks target distribution if called infinitely often

            if (mode == CoercionMode.COERCION_ON) {
                double[] center = parameter.getParameterValues();

//                double balance = 0.0;
                for (int i = 0; i < center.length; ++i) {
                    final double epsilon = scaleFactor * MathUtils.nextGaussian();
                    double x = center[i] + epsilon;
//                    final double old = original[i] - center[i];
//                    balance += epsilon * epsilon - old * old;
                }

//                logHR = 0.5 * balance / (scaleFactor * scaleFactor);
                logHR = 0.0;  // TODO Detailed balance is only met when mode finder run to convergence
            }

            return logHR;

        } else {
            throw new OperatorFailedException("Finished max times");
        }
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return OPERATOR;
    }

    @Override
    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    @Override
    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return scaleFactor;
    }

    public final String getPerformanceSuggestion() {
        return "None";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);

            int maxTimes = xo.getAttribute(MAX_TIMES, Integer.MAX_VALUE);

            MultiDimensionalScalingMM mm = (MultiDimensionalScalingMM) xo.getChild(MMAlgorithm.class);

            return new ModeFindOperator(mm, maxTimes, weight, mode);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a mode finder that always accepts.";
        }

        public Class getReturnType() {
            return ModeFindOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule(MAX_TIMES, true),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(MultiDimensionalScalingMM.class),
        };

    };
}
