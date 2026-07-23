/*
 * NumericalHessianFromGradient.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NumericalHessianFromGradient implements HessianWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradientProvider;

    public NumericalHessianFromGradient(GradientWrtParameterProvider gradientWrtParameterProvider) {
        this.gradientProvider = gradientWrtParameterProvider;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        final int dim = gradientProvider.getDimension();
        double[][] numericalHessian = getNumericalHessianCentral(); //todo: no need to get the full hessian if only need the diagonals
        double[] numericalHessianDiag = new double[dim];

        for (int i = 0; i < dim; i++) {
            numericalHessianDiag[i] = numericalHessian[i][i];
        }
        return numericalHessianDiag;
    }

    @Override
    public double[][] getHessianLogDensity() {
        return getNumericalHessianCentral();
    }

    @Override
    public Likelihood getLikelihood() {
        return gradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return gradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return gradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return gradientProvider.getGradientLogDensity();
    }

    private double[][] getNumericalHessianCentral() {

        final int dim = gradientProvider.getDimension();
        double[][] hessian = new double[dim][dim];

        final double[] oldValues = gradientProvider.getParameter().getParameterValues();

        double[][] gradientPlus = new double[dim][dim];
        double[][] gradientMinus = new double[dim][dim];

        double[] h = new double[dim];
        for (int i = 0; i < dim; i++) {
            h[i] = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(oldValues[i]) + 1.0);
            gradientProvider.getParameter().setParameterValue(i, oldValues[i] + h[i]);
            gradientPlus[i] = gradientProvider.getGradientLogDensity();

            gradientProvider.getParameter().setParameterValue(i, oldValues[i] - h[i]);
            gradientMinus[i] = gradientProvider.getGradientLogDensity();
            gradientProvider.getParameter().setParameterValue(i, oldValues[i]);
        }

        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                hessian[j][i] = hessian[i][j] = (gradientPlus[j][i] - gradientMinus[j][i]) / (4.0 * h[j]) + (gradientPlus[i][j] - gradientMinus[i][j]) / (4.0 * h[i]);
            }
        }

        return hessian;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                GradientWrtParameterProvider.TOLERANCE);
    }
}
