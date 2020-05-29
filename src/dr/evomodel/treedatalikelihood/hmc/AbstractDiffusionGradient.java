/*
 * AbstractDiffusionGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.List;

public abstract class AbstractDiffusionGradient implements GradientWrtParameterProvider, Reportable {

    private final Likelihood likelihood;
    //    private final Parameter parameter;
    private final double lowerBound;
    private final double upperBound;
    protected int offset;

    AbstractDiffusionGradient(Likelihood likelihood, double upperBound, double lowerBound) {
        this.likelihood = likelihood;
//        this.parameter = parameter;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.offset = 0;
    }

    public abstract double[] getGradientLogDensity(double[] gradient);

    public abstract Parameter getRawParameter();

    public void setOffset(int offset) {
        this.offset = offset;
    }

    public abstract ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter();

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

//    @Override
//    public Parameter getParameter() {
//        return parameter;
//    }

    protected Parameter getNumericalParameter() {
        return getParameter();
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                lowerBound, upperBound, TOLERANCE);
    }

//    @Override
//    public String getReport() {
//        return checkNumeric(getGradientLogDensity());
//    }

    String getReportString(double[] analytic, double[] numeric) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric: " + new Vector(numeric) +
                "\n";
    }

    String getReportString(double[] analytic, double[] numeric, double[] numericTrans) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric (no Cholesky): " + new Vector(numeric) +
                "\n" +
                "numeric (with Cholesky): " + new Vector(numericTrans) +
                "\n";
    }

    MultivariateFunction getNumeric() {

        return new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    getNumericalParameter().setParameterValue(i, argument[i]);
                }

                likelihood.makeDirty();
                System.err.println("likelihood in numeric:" + likelihood.getLogLikelihood());
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return lowerBound;
            }

            @Override
            public double getUpperBound(int n) {
                return upperBound;
            }
        };
    }

    String checkNumeric(double[] analytic) {

        System.err.println("Numeric at: \n" + new Vector(getNumericalParameter().getParameterValues()));

        double[] storedValues = getNumericalParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(getNumeric(), storedValues);
        for (int i = 0; i < storedValues.length; ++i) {
            getNumericalParameter().setParameterValue(i, storedValues[i]);
        }

        return getReportString(analytic, testGradient);
    }

    // Default branch specific class
    public static class ParameterDiffusionGradient extends AbstractDiffusionGradient implements Reportable {

        protected final int dim;
        private final BranchSpecificGradient branchSpecificGradient;

        private final Parameter parameter;
        private final Parameter rawParameter;

        ParameterDiffusionGradient(BranchSpecificGradient branchSpecificGradient,
                                   Likelihood likelihood,
                                   Parameter parameter,
                                   Parameter rawParameter,
                                   double upperBound,
                                   double lowerBound) {
            super(likelihood, upperBound, lowerBound);

            this.parameter = parameter;
            this.rawParameter = rawParameter;

            this.branchSpecificGradient = branchSpecificGradient;
            this.dim = parameter.getDimension();

        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        @Override
        public int getDimension() {
            return dim;
        }

        @Override
        public Parameter getRawParameter() {
            return rawParameter;
        }

        @Override
        public double[] getGradientLogDensity() {
            double[] gradient = branchSpecificGradient.getGradientLogDensity();
            return getGradientLogDensity(gradient);
        }

        public double[] getGradientLogDensity(double[] gradient) {
            return extractGradient(gradient);
        }

        private double[] extractGradient(double[] gradient) {
            double[] result = new double[dim];
            for (int i = 0; i < dim; i++) {
                result[i] = gradient[offset + i];
            }
            return result;
        }

        @Override
        public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
            List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derList = branchSpecificGradient.getDerivationParameter();
            assert derList.size() == 1;
            return derList.get(0);
        }

        @Override
        public String getReport() {
            return "Gradient." + rawParameter.getParameterName() + "\n" +
                    super.getReport();
        }

        public static ParameterDiffusionGradient createDriftGradient(BranchSpecificGradient branchSpecificGradient,
                                                                     Likelihood likelihood,
                                                                     Parameter drift) {
            return new ParameterDiffusionGradient(
                    branchSpecificGradient, likelihood,
                    drift, drift,
                    Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        }

        public static ParameterDiffusionGradient createDiagonalAttenuationGradient(BranchSpecificGradient branchSpecificGradient,
                                                                                   Likelihood likelihood,
                                                                                   MatrixParameterInterface attenuation) {
            assert (attenuation instanceof DiagonalMatrix)
                    : "DiagonalAttenuationGradient can only be applied to a DiagonalMatrix.";

            return new ParameterDiffusionGradient(
                    branchSpecificGradient, likelihood,
                    ((DiagonalMatrix) attenuation).getDiagonalParameter(),
                    (DiagonalMatrix) attenuation,
                    Double.POSITIVE_INFINITY, 0.0);
        }
    }
}
