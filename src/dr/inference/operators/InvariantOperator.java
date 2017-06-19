/*
 * InvariantOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.DirtyLikelihoodOperatorParser;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 */
public abstract class InvariantOperator extends SimpleMCMCOperator implements GibbsOperator {

    public InvariantOperator(Parameter parameter, Likelihood likelihood, double weight,
                             boolean checkLikelihood) {
        this.parameter = parameter;
        this.likelihood = likelihood;
        this.checkLikelihood = DEBUG ? true : checkLikelihood;

        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return DirtyLikelihoodOperatorParser.TOUCH_OPERATOR;
    }

    @Override
    public double doOperation() {

        double logLikelihood = 0;
        if (checkLikelihood) {
            if (likelihood != null) {
                logLikelihood = likelihood.getLogLikelihood();
            }
        }

        transform(parameter);

        if (checkLikelihood) {
            if (likelihood != null) {
                double newLogLikelihood = likelihood.getLogLikelihood();

                if (Math.abs(logLikelihood - newLogLikelihood) >  tolerance) {
                    System.err.println("Likelihood is not invariant to transformation:");
                    System.err.println("Before: " + logLikelihood);
                    System.err.println("After : " + newLogLikelihood);
                    System.exit(-1);
                }
            }
        }

        return 0;
    }

    protected abstract void transform(Parameter parameter);

    public int getStepCount() {
        return 1;
    }

    private final Parameter parameter;
    private final Likelihood likelihood;
    private final boolean checkLikelihood;

    private final static boolean DEBUG = false;
    private final static double tolerance = 1E-1;

    public static class Rotation extends InvariantOperator {

        private final boolean translationInvariant;
        private final boolean rotationInvariant;

        public Rotation(Parameter parameter, double weight, Likelihood likelihood,
                        boolean translate, boolean rotate,
                        boolean checkLikelihood) {
            super(parameter, likelihood, weight, checkLikelihood);

            this.translationInvariant = translate;
            this.rotationInvariant = rotate;
        }

        @Override
        protected void transform(Parameter parameter) {

            double[] x = parameter.getParameterValues();

            EllipticalSliceOperator.transformPoint(x, translationInvariant, rotationInvariant, 2);

            final int len = x.length;
            for (int i = 0; i < len; ++i) {
                parameter.setParameterValueQuietly(i, x[i]);
            }
            parameter.fireParameterChangedEvent();
        }
    }
}
