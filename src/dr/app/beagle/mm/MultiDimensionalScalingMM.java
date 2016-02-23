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
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.EllipticalSliceOperator;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

/**
 * Created by msuchard on 12/15/15.
 */
public class MultiDimensionalScalingMM extends MMAlgorithm {

    private final MultiDimensionalScalingLikelihood likelihood;
    private final GaussianProcessRandomGenerator gp;

    private final int P; // Embedding dimension
    private final int Q; // Data dimension

    private double[] XtX = null;
    private double[] D = null;
    private double[] distance = null;

    public MultiDimensionalScalingMM(MultiDimensionalScalingLikelihood likelihood,
                                     GaussianProcessRandomGenerator gp) {
        super();

        this.likelihood = likelihood;
        this.gp = gp;

        this.P = likelihood.getMdsDimension();
        this.Q = likelihood.getLocationCount();

        double[][] precision = gp.getPrecisionMatrix();
        setPrecision(precision);

        double[] mode = null;

        System.err.println("Start: " + printArray(likelihood.getMatrixParameter().getParameterValues()));

        try {
            mode = findMode(likelihood.getMatrixParameter().getParameterValues());
        } catch (NotConvergedException e) {
            e.printStackTrace();
        }

        System.err.println("Final: " + printArray(mode));

        EllipticalSliceOperator.transformPoint(mode, true, true, P);

        System.err.println("Final: " + printArray(mode));

        setParameterValues(likelihood.getMatrixParameter(), mode);
    }

    private void setParameterValues(MatrixParameterInterface mat, double[] values) {

        mat.setAllParameterValuesQuietly(values, 0);
        mat.setParameterValueNotifyChangedAll(0, 0, values[0]); // Fire changed
    }

    private double[] getDistanceMatrix() {
        return likelihood.getObservations();
    }

    private void setPrecision(double[][] matrix) {

        if (!ignoreGP) {

            final int QP = matrix.length;
            if (QP != this.Q * this.P) throw new IllegalArgumentException("Invalid dimensions");

            precision = new double[QP * QP];
            precisionSign = new int[QP * QP];

            precisionStatistics = new double[QP];

            for (int ik = 0; ik < QP; ++ik) {
                double sum = 0.0;
                for (int jl = 0; jl < QP; ++jl) {
                    double value = weightTree * matrix[ik][jl];
                    if (ik != jl) {
                        sum += Math.abs(value);
                    }
                    precisionSign[ik * QP + jl] = (value > 0.0) ? 1 : -1;
                    precision[ik * QP + jl] = Math.abs(value);
                }
                precisionStatistics[ik] = sum;
            }

        }

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
                double norm2 = XtX[i * Q + i] + XtX[j * Q + j] - 2 * XtX[i * Q + j];
                double norm = norm2 > 0.0 ? Math.sqrt(norm2) : 0.0;
                D[j * Q + i] = D[i * Q + j] = Math.max(norm, 1E-3);

                if (Double.isNaN(D[i * Q + j])) {
                    System.err.println("D NaN");
                    System.err.println(XtX[i * Q + i]);
                    System.err.println(XtX[j * Q + j]);
                    System.err.println(2 * XtX[i * Q + j]);
                    System.err.println(norm2);
                    System.err.println(norm);
                    System.exit(-1);
                }
            }
        }

        // Compute update
        for (int i = 0; i < Q; ++i) { // TODO Embarrassingly parallel
            for (int k = 0; k < P; ++k) { // TODO Embarrassingly parallel

                final int ik = i * P + k;
                final int QP = Q * P;

                double numerator = 0.0;
                for (int j = 0; j < Q; ++j) {
                    int add = (i != j)? 1 : 0;
                    double inc = distance[i * Q + j] * (current[i * P + k] - current[j * P + k]) / D[i * Q + j]
                                                + (current[i * P + k] + current[j * P + k]);
                    if (Double.isNaN(inc)) {
                        System.err.println("Bomb at " + i + " " + k + " " + j);
                        System.err.println("Distance = " + distance[i * Q + j]);
                        System.err.println("Ci = " + current[i * P + k]);
                        System.err.println("Cj = " + current[j * P + k]);
                        System.err.println("D = " + D[i * Q + j]);
                        System.exit(-1);
                    }

                    if (precision != null) {
                        for (int l = 0; l < P; ++l) {
                            final int jl = j * P + l;
                            inc += precision[ik * QP + jl] * (current[i * P + k] - precisionSign[ik * QP + jl] * current[j * P + k]);
                        }
                    }

                    numerator += add * inc;
                }
                double denominator = 2 * (Q - 1);

                if (precision != null) {
                    denominator += 2 * precision[ik * QP + ik] + precisionStatistics[ik];
                }

                next[i * P + k] = numerator / denominator;
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

            GaussianProcessRandomGenerator gp =
                    (GaussianProcessRandomGenerator) xo.getChild(GaussianProcessRandomGenerator.class);

            return new MultiDimensionalScalingMM(likelihood, gp);
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
                new ElementRule(GaussianProcessRandomGenerator.class),
        };

        public Class getReturnType() {
            return MultiDimensionalScalingMM.class;
        }
    };

    private double[] precision = null;

    private double[] precisionStatistics = null;
    private int[] precisionSign = null;

    private boolean ignoreGP = false;

    private double weightTree = 1.0;
}
