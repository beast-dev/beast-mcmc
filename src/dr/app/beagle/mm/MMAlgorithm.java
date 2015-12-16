/*
 * MultiDimensionalScalingMM.java
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

package dr.app.beagle.mm;

import dr.app.beagle.multidimensionalscaling.MultiDimensionalScalingLikelihood;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.xml.*;

/**
 * Created by msuchard on 12/15/15.
 */
public class MMAlgorithm {

    public static final double DEFAULT_TOLERANCE = 1E-6;
    public static final int DEFAULT_MAX_ITERATIONS = 1000;

    private final MultiDimensionalScalingLikelihood likelihood;
    private final FullyConjugateMultivariateTraitLikelihood tree;

    public MMAlgorithm(MultiDimensionalScalingLikelihood likelihood,
                       FullyConjugateMultivariateTraitLikelihood tree) {
        this.likelihood = likelihood;
        this.tree = tree;
    }

    public double[] findMode(final double[] startingValue) throws NotConvergedException {
        return findMode(startingValue, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATIONS);
    }

    public double[] findMode(final double[] startingValue, final double tolerance,
                             final int maxIterations) throws NotConvergedException {

        double[] buffer1 = new double[startingValue.length];
        double[] buffer2 = new double[startingValue.length];

        double[] current = buffer1;
        double[] next = buffer2;

        System.arraycopy(startingValue, 0, current, 0, startingValue.length);
        int iteration = 0;

        do {
            mmUpdate(current, next);
            ++iteration;
        } while (convergenceCriterion(next, current) > tolerance && iteration < maxIterations);

        if (iteration >= maxIterations) {
            throw new NotConvergedException();
        }

        return next;
    }

    private void mmUpdate(final double[] current, double[] next) {
        // Do nothing yet
    }

    private double convergenceCriterion(final double[] current, final double[] previous) {
        double norm = 0.0;

        for (int i = 0; i < current.length; ++i) {
            norm += (current[i] - previous[i]) * (current[i] - previous[i]);
        }

        return Math.sqrt(norm);
    }


    class NotConvergedException extends Exception {
        // Nothing interesting
    }

   // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String MDS_STARTING_VALUES = "mdsModeFinder";

        public String getParserName() {
            return MDS_STARTING_VALUES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            return new MMAlgorithm();
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides a mode finder for a MDS model on a tree";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MultiDimensionalScalingLikelihood.class),
                new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
        };

        public Class getReturnType() {
            return MMAlgorithm.class;
        }
    };


}
