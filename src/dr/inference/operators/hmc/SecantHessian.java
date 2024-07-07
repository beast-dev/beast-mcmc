/*
 * SecantHessian.java
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
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
class SecantHessian implements HessianWrtParameterProvider {
    private final int dim;
    private final GradientWrtParameterProvider gradientProvider;
    private final int secantSize;

    private final Secant[] queue;
    private int secantIndex;
    private int secantUpdateCount;
    private double[][] secantHessian;  //TODO: instead store and update the Cholesky decomposition H = LL^T

    SecantHessian(GradientWrtParameterProvider gradientProvider, int secantSize) {
        this.gradientProvider = gradientProvider;
        this.secantSize = secantSize;
        this.dim = gradientProvider.getDimension();

        this.secantHessian = new double[dim][dim]; //TODO: use WrappedVector
        for (int i = 0; i < dim; i++) {
            secantHessian[i][i] = 1.0;
        }

        this.queue = new Secant[secantSize];
        this.secantIndex = 0;
        this.secantUpdateCount = 0;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        double[] diagonalHessian = new double[dim];
        for (int i = 0; i < dim; i++) {
            diagonalHessian[i] = secantHessian[i][i];
        }
        return diagonalHessian;
    }

    @Override
    public double[][] getHessianLogDensity() {
        return secantHessian.clone();
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

    private class Secant { // TODO Inner class because we may change the storage depending on efficiency

        ReadableVector gradient;
        ReadableVector position;

        WrappedVector sk; // sk = dk = x_{k+1} - x_k
        WrappedVector yk; // yk = g_{k+1} - g_k  As defined in 7.18 on page 177 of Numerical Optimization 2nd by Nocedal and Wright
        private double reciprocalInnerProduct;

        Secant(ReadableVector gradient, ReadableVector position) {
            this.gradient = WrappedVector.Utils.copy(gradient);
            this.position = WrappedVector.Utils.copy(position);

            this.sk = new WrappedVector.Raw(new double[position.getDim()]);
            this.yk = new WrappedVector.Raw(new double[position.getDim()]);
            this.reciprocalInnerProduct = 0.0;
        }

        double getPosition(int index) { return position.get(index); }

        double getGradient(int index) {
            return gradient.get(index);
        }

        void updateSkYk(Secant newSecant) {

            reciprocalInnerProduct = 0.0;

            for (int i = 0; i < dim; i++) {
                double ski = newSecant.getPosition(i) - getPosition(i);
                double yki = newSecant.getGradient(i) - getGradient(i);
                sk.set(i, ski);
                yk.set(i, yki);
                reciprocalInnerProduct += sk.get(i) * yk.get(i);
            }

            reciprocalInnerProduct = 1.0 / reciprocalInnerProduct;
        }

        double getReciprocalInnerProduct() {
            return reciprocalInnerProduct;
        }

        double getYk(int i) {
            return yk.get(i);
        }

        double getSk(int i) {
            return sk.get(i);
        }
    }

    public void storeSecant(ReadableVector gradient, ReadableVector position) {
        queue[secantIndex] = new Secant(gradient, position);

        if (secantUpdateCount == 0) {
            initializeSecantHessian(queue[secantIndex]);
        } else {
            final int lastSecantIndex = (secantIndex + secantSize - 1) % queue.length;
            updateSecantHessian(queue[secantIndex], queue[lastSecantIndex]);
        }

        secantIndex = (secantIndex + 1) % queue.length;
        ++secantUpdateCount;
    }

    private void initializeSecantHessian(Secant secant) {
            /*
            Equation 6.20 on page 143 of Numerical Optimization 2nd by Nocedal and Wright
             */
        double numerator = 0.0;
        double denominator = 0.0;
        for (int i = 0; i < dim; i++) {
            numerator += secant.getGradient(i) * secant.getPosition(i);
            denominator += secant.getGradient(i) * secant.getGradient(i);
        }

        final double fraction = (denominator == 0 ? 1.0 : numerator / denominator);

        for (int i = 0; i < dim; i++) {
            secantHessian[i][i] = fraction;
        }
    }

    private void updateSecantHessian(Secant lastSecant, Secant currentSecant) {
            /*
            Implementation of formula 11.17 on Page 281 of Ken Lange's Optimization book.
             */
//        double[] d = new double[dim];
//        double[] g = new double[dim];
        double[] Hd = new double[dim];
//        double b = 0.0;
        double c = 0.0;

//        for (int i = 0; i < dim; i++) {
//            d[i] = currentSecant.getPosition(i) - lastSecant.getPosition(i);
//            g[i] = currentSecant.getGradient(i) - lastSecant.getGradient(i);
//            b += d[i] * g[i];
//        }
//        b = 1.0 / b;

        lastSecant.updateSkYk(currentSecant);
        double b = lastSecant.getReciprocalInnerProduct();

        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j =0; j < dim; j++) {
//                sum += secantHessian[i][j] * d[j];
                sum += secantHessian[i][j] * lastSecant.getSk(j);
            }
            Hd[i] = sum;
//            c += d[i] * Hd[i];
            c += lastSecant.getSk(i) * Hd[i];
        }
        c = -1.0 / c;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
//                secantHessian[i][j] += b * g[i] * g[j] + c * Hd[i] * Hd[j];
                secantHessian[i][j] += b * lastSecant.getYk(i) * lastSecant.getYk(j) + c * Hd[i] * Hd[j];
            }
        }
    }
}
