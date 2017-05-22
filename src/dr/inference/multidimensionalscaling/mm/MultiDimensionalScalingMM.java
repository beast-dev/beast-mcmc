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

package dr.inference.multidimensionalscaling.mm;

import dr.inference.multidimensionalscaling.MultiDimensionalScalingLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.EllipticalSliceOperator;
import dr.math.distributions.GaussianProcessRandomGenerator;
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

    final private double tolerance;

    public MultiDimensionalScalingLikelihood getLikelihood() { return likelihood; }

    public GaussianProcessRandomGenerator getGaussianProcess() { return gp; }

    public double getTolerance() { return tolerance; }

    public MultiDimensionalScalingMM(MultiDimensionalScalingLikelihood likelihood,
                                     GaussianProcessRandomGenerator gp,
                                     double tolerance) {
        super();

        this.likelihood = likelihood;
        this.gp = gp;

        this.P = likelihood.getMdsDimension();
        this.Q = likelihood.getLocationCount();

        this.tolerance = tolerance;
    }

    public void run() {
        run(100000);
    }

    public void run(final int maxIterations) {

        if (maxIterations == 0) return;

        if (gp != null) {
            double[][] precision = gp.getPrecisionMatrix();
            setPrecision(precision);
        }

//        System.err.println("");

//        System.err.println("weight: " + 1.0 / likelihood.getMDSPrecision());
        weightTree =  1.0 / likelihood.getMDSPrecision();

        double[] start = likelihood.getMatrixParameter().getParameterValues();

        System.err.println("Start: " + printArray(start));
        double penaltyStart = printLogObjective();

//        EllipticalSliceOperator.transformPoint(mode, true, true, P);
//        setParameterValues(likelihood.getMatrixParameter(), mode);
//        printLogObjective();

        double[] mode = null;

        try {
            mode = findMode(likelihood.getMatrixParameter().getParameterValues(),
                    tolerance, maxIterations);

        } catch (NotConvergedException e) {
            e.printStackTrace();
        }

//        System.err.println("Final: " + printArray(mode));

        setParameterValues(likelihood.getMatrixParameter(), mode);
        double penaltyEnd = printLogObjective();

        System.err.println("Move: " + penaltyStart + " -> " + penaltyEnd + " : " + (penaltyEnd - penaltyStart));
//
//        if (penaltyStart - penaltyEnd > 1E-1) {

//        if (penaltyEnd < penaltyStart) {
//            throw new RuntimeException("MM moved up-hill\n\tStart: " + penaltyStart + "\n\tEnd  : " + penaltyEnd);
//            System.err.println("Revert:  MM moved up-hill?");
//            setParameterValues(likelihood.getMatrixParameter(), start);
//            double penaltyRevert = printLogObjective();
//            System.err.println("End: " + penaltyEnd);
//            System.err.println("revert : " + penaltyRevert);
           // throw new RuntimeException("out");
//
//        }

//        EllipticalSliceOperator.transformPoint(mode, true, true, P);
//
//        System.err.println("Final: " + printArray(mode));
//
//        setParameterValues(likelihood.getMatrixParameter(), mode);
//        printLogObjective();

//        throw new RuntimeException("done");
    }

    private double printLogObjective() {
        double logLike = likelihood.getLogLikelihood();
        double logPenalty = gp.getLikelihood().getLogLikelihood();
        double total = logLike;
        if (weightTree != 0.0) {
            total += logPenalty;
        }
        System.err.println("obj: " + total + " = " + logLike + " + " + logPenalty);
//        return logPenalty;
        return total;
    }

    private void setParameterValues(MatrixParameterInterface mat, double[] values) {

//        for (int i = 0; i < values.length; ++i) {
//            mat.setValue(i, values[i]);
//        }
        mat.setAllParameterValuesQuietly(values, 0);
        mat.setParameterValueNotifyChangedAll(0, 0, values[0]); // Fire changed
//        for (int i = 0; i < mat.getUniqueParameterCount(); ++i) {
//            mat.getUniqueParameter(0).fireParameterChangedEvent();
//        }
//        mat.getUniqueParameter(0);
    }

    private double[] getDistanceMatrix() {
        return likelihood.getObservations();
    }

    private void setPrecision(double[][] matrix) {

        if (!ignoreGP) {

            final int QP = matrix.length;
            if (QP != this.Q * this.P) throw new IllegalArgumentException("Invalid dimensions");

            precision = matrix;
//            precision = new double[QP * QP];
//            precisionSign = new int[QP * QP];

            precisionStatistics = new double[QP];

            for (int ik = 0; ik < QP; ++ik) {
                double sum = 0.0;
                for (int jl = 0; jl < QP; ++jl) {
//                    double value = weightTree * matrix[ik][jl];
                    if (ik != jl) {
                        sum += Math.abs(precision[ik][jl]);
                    }
//                    precisionSign[ik * QP + jl] = (value > 0.0) ? 1 : -1;
//                    precision[ik * QP + jl] = Math.abs(value);
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
                D[j * Q + i] = D[i * Q + j] = Math.max(norm, 1E-10);

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
//                final int QP = Q * P;

                double numerator = 0.0;
                for (int j = 0; j < Q; ++j) {
                    double inc = 0.0;
                    if (i != j) {
//                        int add = (i != j) ? 1 : 0;
//                        double inc = add * distance[i * Q + j] * (current[i * P + k] - current[j * P + k]) / D[i * Q + j]
//                                + (current[i * P + k] + current[j * P + k]);
                        inc += distance[i * Q + j] * (current[i * P + k] - current[j * P + k]) / D[i * Q + j]
                                                                        + (current[i * P + k] + current[j * P + k]);
                    }
//                    inc = 0.0; // TODO Remove!

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
                            final double prec = precision[ik][jl];
                            final int sign = (prec > 0.0) ? 1 : -1;
                            inc += weightTree * Math.abs(prec) * (current[i * P + k] - sign * current[j * P + l]);
                        }
                    }

                    numerator += inc;
                }
                double denominator = 2 * (Q - 1);

//                denominator = 0.0; // TODO Remove

                if (precision != null) {
                    denominator += weightTree * (2 * precision[ik][ik] + precisionStatistics[ik]);
                }

                next[i * P + k] = numerator / denominator;
            }
        }

        // Force translation, rotation, reflection symmetry
        EllipticalSliceOperator.transformPoint(next, true, true, P);
    }

   // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String MDS_STARTING_VALUES = "mdsModeFinder";
        public static final String TOLERANCE = "tolerance";

        public String getParserName() {
            return MDS_STARTING_VALUES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MultiDimensionalScalingLikelihood likelihood =
                    (MultiDimensionalScalingLikelihood) xo.getChild(MultiDimensionalScalingLikelihood.class);

            GaussianProcessRandomGenerator gp =
                    (GaussianProcessRandomGenerator) xo.getChild(GaussianProcessRandomGenerator.class);

            double tolerance = xo.getAttribute(TOLERANCE, 1E-3);

            MultiDimensionalScalingMM mm = new MultiDimensionalScalingMM(likelihood, gp, tolerance);
            mm.run();
            return mm;
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
                new ElementRule(GaussianProcessRandomGenerator.class, true),
                AttributeRule.newDoubleRule(TOLERANCE, true),
        };

        public Class getReturnType() {
            return MMAlgorithm.class;
        }
    };

    private double[][] precision = null;

    private double[] precisionStatistics = null;
//    private int[] precisionSign = null;

    private boolean ignoreGP = false;

    private double weightTree;  // TODO Formally compute
}
