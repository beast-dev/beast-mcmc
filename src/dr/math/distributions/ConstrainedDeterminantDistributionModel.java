/*
 * ConstrainedDeterminantDistributionModel.java
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

package dr.math.distributions;


import dr.inference.model.GradientProvider;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class ConstrainedDeterminantDistributionModel implements MultivariateDistribution, GradientProvider {

    private final int dim;
    private final double shape;

    public ConstrainedDeterminantDistributionModel(double shape, int dim) {
        this.dim = dim;
        this.shape = shape;
    }


    @Override
    public int getDimension() {
        return dim * dim;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    private double[] gradLogPdf(double[] x) {
        DenseMatrix64F X = DenseMatrix64F.wrap(dim, dim, x);

        DenseMatrix64F Xinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(X, Xinv);

        CommonOps.scale(shape, Xinv);
        CommonOps.transpose(Xinv);
        return Xinv.getData();
    }

    @Override
    public double logPdf(double[] x) {
        DenseMatrix64F X = DenseMatrix64F.wrap(dim, dim, x);
        double det = CommonOps.det(X);

        return shape * Math.log(Math.abs(det)); //TODO: normalizing constant
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("not implemented");
    }

    @Override
    public String getType() {
        return "ConstrainedDeterminant";
    }
}
