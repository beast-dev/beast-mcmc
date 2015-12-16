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
public class MultiDimensionalScalingMM extends MMAlgorithm {

    private final MultiDimensionalScalingLikelihood likelihood;
    private final FullyConjugateMultivariateTraitLikelihood tree;

    private final int P; // Embedding dimension
    private final int Q; // Data dimension

    private double[] XtX = null;
    private double[] D = null;
    private double[] distance = null;

    public MultiDimensionalScalingMM(MultiDimensionalScalingLikelihood likelihood, 
                                     FullyConjugateMultivariateTraitLikelihood tree) {
        super();
        this.likelihood = likelihood;
        this.tree = tree;

        this.P = likelihood.getMdsDimension();
        this.Q = likelihood.getLocationCount();

    }

    private double[] getDistanceMatrix() {
        double[] distance = new double[Q * Q];
        for (int i = 0; i < Q; ++i) {
            for (int j = 0; j < Q; ++j) {
                distance[i * Q + j] = 0.0; // TODO
            }
        }
        return distance;
    }

    protected void mmUpdate(final double[] current, double[] next) {

        if (XtX == null) {
            XtX = new double[Q * Q];
        }

        if (D == null) {
            D = new double[Q * Q];
            for (int i = 0; i < Q; ++i) {
                D[i * Q + i] = 1.0; // To protect against divide-by-zero
            }
        }

        if (distance == null) {
            distance = getDistanceMatrix();
        }

        // Compute XtX
        for (int i = 0; i < Q; ++i) {
            for (int j = i; j < Q; ++j) {
                double innerProduct = 0.0;
                for (int k = 0; k < P; ++k) {
                    innerProduct += current[i * P + k] * current[j * P + k];
                }
                XtX[j * Q + i] = XtX[i * Q + j] = innerProduct;
            }
        }

        // Compute D
        for (int i = 0; i < Q; ++i) {
            for (int j = i + 1; j < Q; ++j) { // TODO XtX is not a necessary intermediate
                D[j * Q + i] = D[i * Q + j] = Math.sqrt(XtX[i * Q + i] + XtX[j * Q + j] - 2 * XtX[i * Q + j]);
            }
        }

        // Compute update
        for (int i = 0; i < Q; ++i) { // TODO Embarrassingly parallel
            for (int k = 0; k < P; ++k) { // TODO Embarrassingly parallel
                double sum = 0.0;
                for (int j = 0; j < Q; ++j) {
                    int add = (i != j)? 1 : 0;
                    sum += add * (
                            distance[i * Q + j] * (current[i * P + k] - current[j * P + k]) / D[i * Q + j]
                            + (current[i * P + k] + current[j * P + k])

                    );
                }
                next[i * P + k] = sum / (2 * (Q - 1));
            }
        }

        // Force translation, rotation, reflection symmetry
        for (int k = 0; k < P; ++k) {
            current[0 * P + k] = 0.0; // translation
            if (k < P - 1) {
                current[1 * P + k] = 0.0; // rotation
            } else {
                current[1 * P + k] = Math.abs(current[1 * P + k]); // reflection
            }
        }
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

            MultiDimensionalScalingLikelihood likelihood =
                    (MultiDimensionalScalingLikelihood) xo.getChild(MultiDimensionalScalingLikelihood.class);
            FullyConjugateMultivariateTraitLikelihood tree =
                    (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

            return new MultiDimensionalScalingMM(likelihood, tree);
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
            return MultiDimensionalScalingMM.class;
        }
    };
}
