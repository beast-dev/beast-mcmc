/*
 * TreeTraitLogger.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.model;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;

import static java.lang.Math.sqrt;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public abstract class AbstractVarianceProportionStatistic extends Statistic.Abstract implements Reportable {

    protected final TreeModel tree;
    protected final RepeatedMeasuresTraitDataModel dataModel;
    protected final TreeDataLikelihood treeLikelihood;

    private DenseMatrix64F diffusionProportion;

    protected DenseMatrix64F diffusionComponent;
    protected DenseMatrix64F samplingComponent;

    private final MatrixRatios ratio;

    protected final int dimTrait;

    public AbstractVarianceProportionStatistic(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                               RepeatedMeasuresTraitDataModel dataModel,
                                               MatrixRatios ratio) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.dataModel = dataModel;
        this.dimTrait = dataModel.getTraitDimension();
        this.diffusionProportion = new DenseMatrix64F(dimTrait, dimTrait);
        this.diffusionComponent = new DenseMatrix64F(dimTrait, dimTrait);
        this.samplingComponent = new DenseMatrix64F(dimTrait, dimTrait);

        this.ratio = ratio;
    }

    @Override
    public String getReport() {
        Matrix mat = new Matrix(dimTrait, dimTrait);
        for (int i = 0; i < dimTrait; i++) {
            int offset = dimTrait * i;
            for (int j = 0; j < dimTrait; j++) {
                mat.set(i, j, getStatisticValue(offset + j));
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Variance proportion statistic: " + ratio.name());
        sb.append("\n");
        sb.append("stat value = ");
        sb.append(mat);
        sb.append("\n\n");
        return sb.toString();
    }

    public enum MatrixRatios {
        ELEMENT_WISE {
            @Override
            void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix,
                                DenseMatrix64F destination) {
                int dim = destination.numRows;

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {

                        double n = Math.abs(numeratorMatrix.get(i, j));
                        double d = Math.abs(otherMatrix.get(i, j));

                        if (n == 0 && d == 0) {
                            destination.set(i, j, 0);
                        } else {
                            destination.set(i, j, n / (n + d));
                        }

                    }
                }

            }
        },
        SYMMETRIC_DIVISION {
            @Override
            void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix, DenseMatrix64F destination)
                    throws IllegalDimension {


                //TODO: implement for eigendecomposition with DensMatrix64F

                throw new RuntimeException(SYMMETRIC_DIVISION + " not yet implemented.");


//                int dim = destination.numRows;

//
//                Matrix M1 = numeratorMatrix.add(otherMatrix); //M1 = numeratorMatrix + otherMatrix
//                Matrix M2 = getMatrixSqrt(M1, true); //M2 = inv(sqrt(numeratorMatrix + otherMatrix))
//                Matrix M3 = M2.product(numeratorMatrix.product(M2));//M3 = inv(sqrt(numeratorMatrix + otherMatrix)) *
//                //                                            numeratorMatrix * inv(sqrt(numeratorMatrix + otherMatrix))
//                for (int i = 0; i < dim; i++) {
//                    for (int j = 0; j < dim; j++) {
//                        destination.set(i, j, M3.component(i, j));
//                    }
//                }

            }
        },
        CO_HERITABILITY {
            @Override
            void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix,
                                DenseMatrix64F destination) {

                for (int i = 0; i < destination.numRows; i++) {

                    double val = numeratorMatrix.get(i, i) / (numeratorMatrix.get(i, i) + otherMatrix.get(i, i));
                    destination.set(i, i, val);
                    for (int j = i + 1; j < destination.numRows; j++) {

                        double rg = numeratorMatrix.get(i, j);
                        double vi = numeratorMatrix.get(i, i) + otherMatrix.get(i, i);
                        double vj = numeratorMatrix.get(j, j) + otherMatrix.get(j, j);

                        val = rg / sqrt(vi * vj);

                        destination.set(i, j, val);
                        destination.set(j, i, val);

                    }
                }
            }
        };

        abstract void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix,
                                     DenseMatrix64F destination)
                throws IllegalDimension;
    }

    private void updateDiffsionProportion() throws IllegalDimension {
        updateVarianceComponents();
        ratio.setMatrixRatio(diffusionComponent, samplingComponent, diffusionProportion);
    }

    abstract protected void updateVarianceComponents();

    //TODO: Move method below to a different class
    //TODO: implement this for DenseMatrix64F rather than Matrix
    private static Matrix getMatrixSqrt(Matrix M, Boolean invert) {
        DoubleMatrix2D S = new DenseDoubleMatrix2D(M.toComponents());
        RobustEigenDecomposition eigenDecomp = new RobustEigenDecomposition(S, 100);
        DoubleMatrix1D eigenValues = eigenDecomp.getRealEigenvalues();
        int dim = eigenValues.size();
        for (int i = 0; i < dim; i++) {
            double value = sqrt(eigenValues.get(i));
            if (invert) {
                value = 1 / value;
            }
            eigenValues.set(i, value);
        }

        DoubleMatrix2D eigenVectors = eigenDecomp.getV();
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                eigenVectors.set(i, j, eigenVectors.get(i, j) * eigenValues.get(j));

            }
        }
        DoubleMatrix2D storageMatrix = new DenseDoubleMatrix2D(dim, dim);
        eigenVectors.zMult(eigenDecomp.getV(), storageMatrix, 1, 0, false, true);


        return new Matrix(storageMatrix.toArray());

    }


    @Override
    public int getDimension() {
        return dimTrait * dimTrait;
    }

    @Override
    public String getDimensionName(int dim) {
        int row = dim / dimTrait;
        int col = dim - row * dimTrait;
        return getStatisticName() + (row + 1) + (col + 1);
    }

    @Override
    public double getStatisticValue(int dim) {

        boolean needToUpdate = needToUpdate(dim);

        if (needToUpdate) {

            try {
                updateDiffsionProportion();
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

        }


        int d1 = dim / dimTrait;
        int d2 = dim - d1 * dimTrait;
        return diffusionProportion.get(d1, d2);

    }

    abstract protected boolean needToUpdate(int dim);
}
