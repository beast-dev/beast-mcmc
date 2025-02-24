/*
 * GeodesicHamiltonianMonteCarloOperator.java
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
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inferencexml.operators.hmc.GeodesicHamiltonianMonteCarloOperatorParser;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;


public class GeodesicHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    public GeodesicHamiltonianMonteCarloOperator(AdaptationMode mode, double weight,
                                                 GradientWrtParameterProvider gradientProvider, Parameter parameter,
                                                 Transform transform, Parameter maskParameter, Options runtimeOptions,
                                                 MassPreconditioner preconditioner,
                                                 ManifoldProvider manifolds) {
        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner);
        this.leapFrogEngine = new GeodesicLeapFrogEngine(parameter,
                getDefaultInstabilityHandler(), preconditioning, mask, manifolds);
    }

    public void addManifolds(ManifoldProvider manifolds) {
        ((GeodesicLeapFrogEngine) leapFrogEngine).addManifolds(manifolds);
    }

    @Override
    public String getOperatorName() {
        return "GeodesicHMC(" + parameter.getParameterName() + ")";
    }

    @Override
    public String getReport() {

        MatrixParameterInterface matParam = (MatrixParameterInterface) parameter;
        int k = matParam.getColumnDimension();
        int p = matParam.getRowDimension();

        StringBuilder sb = new StringBuilder("operator: " + GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME);
        sb.append("\n");
        sb.append("\toriginal position:\n");
        Matrix originalPosition = new Matrix(matParam.getParameterAsMatrix());
        sb.append(originalPosition.toString(2));

        double[] momentum = new double[parameter.getDimension()];

        for (int i = 0; i < momentum.length; i++) { //Need some deterministic way to assign momentum for test
            momentum[i] = i;
        }

        Matrix originalMomentum = new Matrix(p, k);
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < k; j++) {
                originalMomentum.set(i, j, momentum[i + j * p]);
            }
        }
        sb.append("\toriginal momentum (unprojected):\n");
        sb.append(originalMomentum.toString(2));

        WrappedVector wrappedMomentum = new WrappedVector.Raw(momentum);
        double hastings;
        try {
            hastings = leapFrogGivenMomentum(wrappedMomentum);
        } catch (NumericInstabilityException e) {
            e.printStackTrace();
            throw new RuntimeException("HMC failed");
        }

        Matrix finalPosition = new Matrix(matParam.getParameterAsMatrix());
        sb.append("\n");
        sb.append("\tfinal position:\n");
        sb.append(finalPosition.toString(2));
        sb.append("\n");
        sb.append("\thastings ratio: " + hastings + "\n\n");


        return sb.toString();
    }

//    public void setOrthogonalityStructure(ArrayList<ArrayList<Integer>> oldOrthogonalityStructure) {
//        ((GeodesicLeapFrogEngine) leapFrogEngine).setOrthogonalityStructure(oldOrthogonalityStructure);
//    }


    public static class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

        //        private final DenseMatrix64F positionMatrix;
//        private final DenseMatrix64F innerProduct;
//        private final DenseMatrix64F innerProduct2;
//        private final DenseMatrix64F projection;
//        private final DenseMatrix64F momentumMatrix;
//        private final int nRows;
//        private final int nCols;

        //        private final int[] subRows;
//        private final int[] subColumns;


        private ManifoldProvider manifoldProvider;

        GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                               MassPreconditioner preconditioning, double[] mask,
                               ManifoldProvider manifoldProvider) {
            super(parameter, instabilityHandler, preconditioning, mask);
            this.manifoldProvider = manifoldProvider;


//            this.subRows = parseSubRowsFromMask();
//            this.subColumns = parseSubColumnsFromMask();


//            orthogonalityStructure.add(subColumns);

//            this.nRows = subRows.length;

//            if (mask != null) checkMask(subRows, subColumns);
//            this.nCols = subColumns.length;
//            this.positionMatrix = new DenseMatrix64F(nCols, nRows);
//            this.innerProduct = new DenseMatrix64F(nCols, nCols);
//            this.innerProduct2 = new DenseMatrix64F(nCols, nCols);
//            this.projection = new DenseMatrix64F(nCols, nRows);
//            this.momentumMatrix = new DenseMatrix64F(nCols, nRows);
        }


//        private int[] parseSubColumnsFromMask() {
//
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            ArrayList<Integer> subArray = new ArrayList<Integer>();
//
//            for (int col = 0; col < originalColumns; col++) {
//                int offset = col * originalRows;
//                for (int row = 0; row < originalRows; row++) {
//                    int ind = offset + row;
//                    if (mask == null || mask[ind] == 1.0) {
//                        subArray.add(col);
//                        break;
//                    }
//                }
//            }
//
//            int[] subColumns = new int[subArray.size()];
//            for (int i = 0; i < subColumns.length; i++) {
//                subColumns[i] = subArray.get(i);
//            }
//
//            return subColumns;
//        }
//
//        private int[] parseSubRowsFromMask() {
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            ArrayList<Integer> subArray = new ArrayList<Integer>();
//
//            for (int row = 0; row < originalRows; row++) {
//                for (int col = 0; col < originalColumns; col++) {
//                    int ind = col * originalRows + row;
//                    if (mask == null || mask[ind] == 1.0) {
//                        subArray.add(row);
//                        break;
//                    }
//                }
//            }
//
//            int[] subRows = new int[subArray.size()];
//            for (int i = 0; i < subRows.length; i++) {
//                subRows[i] = subArray.get(i);
//            }
//
//            return subRows;
//        }

//        private void checkMask(int[] rows, int[] cols) {
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            int subRowInd = 0;
//            int subColInd = 0;
//
//            Boolean isSubRow;
//            Boolean isSubCol;
//
//            for (int row = 0; row < originalRows; row++) {
//                if (subRowInd < rows.length && row == rows[subRowInd]) {
//                    isSubRow = true;
//                    subRowInd++;
//                } else {
//                    isSubRow = false;
//                }
//
//                subColInd = 0;
//
//                for (int col = 0; col < originalColumns; col++) {
//                    if (subColInd < cols.length && col == cols[subColInd]) {
//                        isSubCol = true;
//                        subColInd++;
//                    } else {
//                        isSubCol = false;
//                    }
//
//                    int ind = originalRows * col + row;
//
//                    if (isSubCol && isSubRow) {
//                        if (mask[ind] != 1.0) {
//                            throw new RuntimeException("mask is incompatible with " +
//                                    GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME +
//                                    ". All elements in sub-matrix must be set to 1.");
//                        }
//                    } else {
//                        if (mask[ind] != 0.0) {
//                            throw new RuntimeException("mask is incompatible with " +
//                                    GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME +
//                                    ". All elements outside of sub-matrix must be set to 0.");
//                        }
//                    }
//
//                }
//            }
//        }

        public void addManifolds(ManifoldProvider manifoldProvider) {
            this.manifoldProvider = manifoldProvider;
        }


        @Override
        public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                   double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {
            super.updateMomentum(position, momentum, gradient, functionalStepSize);
            projectMomentum(momentum, position);

        }

        @Override
        public void updatePosition(double[] position, WrappedVector momentum,
                                   double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

            manifoldProvider.updatePositionAndMomentum(position, momentum, functionalStepSize);
            parameter.setAllParameterValuesQuietly(position);
            parameter.fireParameterChangedEvent();
        }

        @Override
        public void projectMomentum(double[] momentum, double[] position) {
            manifoldProvider.projectTangent(momentum, position);
        }
    }
}
