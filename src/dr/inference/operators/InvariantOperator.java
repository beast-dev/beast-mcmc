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

/**
 * @author Marc Suchard
 */
public abstract class InvariantOperator extends SimpleMCMCOperator implements GibbsOperator {

    private InvariantOperator(Parameter parameter, Likelihood likelihood, double weight,
                             boolean checkLikelihood) {
        this.parameter = parameter;
        this.likelihood = likelihood;
        this.checkLikelihood = checkLikelihood || DEBUG;

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

        if (pathParameter == 1.0) {
            transform(parameter);
        }

        if (checkLikelihood) {
            if (likelihood != null) {
                double newLogLikelihood = likelihood.getLogLikelihood();

                if (Math.abs(logLikelihood - newLogLikelihood) >  tolerance) {
                    String sb = "Likelihood is not invariant to transformation:\n" +
                            "Before: " + logLikelihood + "\n" +
                            "After : " + newLogLikelihood + "\n";
                    throw new RuntimeException(sb);
                }
            }
        }

        return 0;
    }

    @Override
    public void setPathParameter(double beta) {
        pathParameter = beta;
    }

    private double pathParameter = 1.0;

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
        private final int dim;

        public Rotation(Parameter parameter, int dim,
                        double weight, Likelihood likelihood,
                        boolean translate, boolean rotate,
                        boolean checkLikelihood) {
            super(parameter, likelihood, weight, checkLikelihood);

            this.dim = dim;
            this.translationInvariant = translate;
            this.rotationInvariant = rotate;
        }

        @Override
        protected void transform(Parameter parameter) {

            double[] x = parameter.getParameterValues();

            EllipticalSliceOperator.transformPoint(x, translationInvariant, rotationInvariant, dim);

            final int len = x.length;
            for (int i = 0; i < len; ++i) {
                parameter.setParameterValueQuietly(i, x[i]);
            }
            parameter.fireParameterChangedEvent();
        }
    }
}
